#define NOMINMAX
#include <vector>
#include <stdio.h>
#include <string.h>
#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <sstream>

#include <proj.h>

#include <colmap/base/reconstruction.h>
#include <colmap/feature/extraction.h>
#include <colmap/feature/matching.h>

#include "colmap/base/reconstruction.h"
#include "colmap/controllers/automatic_reconstruction.h"
#include "colmap/controllers/bundle_adjustment.h"
#include "colmap/controllers/hierarchical_mapper.h"
#include "colmap/exe/gui.h"
#include "colmap/exe/sfm.h"
#include "colmap/util/bitmap.h"
#include "colmap/exe/model.h"
#include "colmap/util/misc.h"
#include "colmap/util/opengl_utils.h"
#include "colmap/util/option_manager.h"
#include "colmap/base/gps.h"
#include "colmap/util/ply.h"

#include "colmap/base/undistortion.h"
#include "mvs/fusion.h"
#include "mvs/meshing.h"
#include "mvs/patch_match.h"
#include "util/misc.h"

#include "boost/format.hpp"
#include "exif.h"

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include "exiv2/exiv2.hpp"

/*
 * Todo
 * 1. Confirm the pose (both position and orientation)
 * 2. Fix the intrinsic matrix
 * 3. Export ply pointclouds
 *
 * Advanced
 * 1. "Dense" sfm
 */

#define ACCEPT_USR_OF_DEPRECATED_PROJ_API_H 1
namespace fs = boost::filesystem;
using namespace colmap;

 // parse the gps info of [photoPath] into xyz
Eigen::Vector3d getImageLocation_revised(const std::string& v_path, PJ* v_transformer, const Eigen::Vector3d& v_origin)
{
	easyexif::EXIFInfo result;
	{
		FILE* fp = fopen(v_path.c_str(), "rb");
		if (!fp) {
			printf("Can't open file.\n");
			throw;
		}
		fseek(fp, 0, SEEK_END);
		unsigned long fsize = ftell(fp);
		rewind(fp);
		unsigned char* buf = new unsigned char[fsize];
		if (fread(buf, 1, fsize, fp) != fsize) {
			printf("Can't read file.\n");
			delete[] buf;
			throw;
		}
		fclose(fp);

		int code = result.parseFrom(buf, fsize);
		delete[] buf;
		if (code) {
			printf("Error parsing EXIF: code %d\n", code);
			throw;
		}
	}

	double longtitude = result.GeoLocation.Longitude;
	double latitude = result.GeoLocation.Latitude;
	double altitude = result.GeoLocation.Altitude;

	std::cout << latitude << ' ' << longtitude << ' ' << altitude << std::endl;

	PJ_COORD target_ = proj_trans(v_transformer, PJ_FWD, proj_coord(longtitude, latitude, 0, 0));
	Eigen::Vector3d target(
		target_.lp.lam,
		target_.lp.phi,
		altitude
	);
	target -= v_origin;
	std::cout << target[0] << ' ' << target[1] << ' ' << target[2] << std::endl;
	return target;
}

Eigen::Matrix<double, 4, 4> ReadMatrixFromFile(const std::string& filename) {
	Eigen::Matrix<double, 4, 4> matrix = Eigen::Matrix<double, 4, 4>::Identity();
	std::ifstream file(filename);

	if (!file.is_open()) {
		std::cerr << "cannot open file: " << filename << std::endl;
		throw std::exception();
	}

	for (int i = 0; i < 4; ++i) {
		std::string line;
		if (std::getline(file, line)) {
			std::istringstream line_stream(line);
			for (int j = 0; j < 4; ++j) {
				double value;
				if (line_stream >> value) {
					matrix(i, j) = value;
				}
			}
		}
	}
	return matrix;
}

void checkFolder(const fs::path& folder) {
	if (fs::is_directory(folder)) {
		fs::remove_all(folder);
	}
	fs::create_directories(folder);
}

void safeCheckFolder(const fs::path& folder) {
	if (!fs::is_directory(folder)) {
		fs::create_directories(folder);
	}
}

void ExportPLYwithNormal(Reconstruction& reconstruction, const std::string& path) {
	const auto ply_points = reconstruction.ConvertToPLY();
	const bool kWriteNormal = true;
	const bool kWriteRGB = false;
	WriteBinaryPlyPoints(path, ply_points, kWriteNormal, kWriteRGB);
}

int runcolmap_with_align()
{
	//fs::path output_root = "H:/data/L7_part";
	//const std::string db_path = "H:/data/L7_part/database.db";  // Canbe blank
	//const std::string img_path = "H:/data/L7_part/images";      // Path of the input image

	fs::path output_root = "H:/Data/runcolmapData";
	const std::string db_path = "H:/Data/runcolmapData/database.db";  // Canbe blank
	const std::string img_path = "H:/Data/runcolmapData/images";      // Path of the input image

	// global config
	bool isRunningAlignment = true;

	// Some config option
	// Pro and Prior option
	bool use_4547 = true;
	// const Eigen::Vector3d origin_point(493260.00, 2492700.00, 0);
	// const double focal_length_in_pixel = 9004;

	const Eigen::Vector3d origin_point(492940.0287, 2493582.533, 45.45);
	const double focal_length_in_pixel = 8192;

	// Recon option
	bool isHighQuality = false;
	bool isMultipleModels = false;
	bool isUseGpu = true;
	bool isSingleCamera = false;
	bool isGreaterBa = true;
	std::int16_t max_image_size = 8192;

	// Load exif and rewrite the image to rectify the orientation
	LOG(INFO) << "===== Load exif and rewrite the image to rectify the orientation =====";
	if (!fs::exists(output_root / "imgs"))
	{
		fs::directory_iterator it_dir(img_path), it_dir_end;

		std::vector<std::string> img_paths;

		for (; it_dir != it_dir_end; ++it_dir)
			if (it_dir->path().extension() == ".png" || it_dir->path().extension() == ".PNG" || it_dir->path().extension() == ".jpg" || it_dir->path().extension() == ".JPG")
				img_paths.emplace_back(it_dir->path().filename().string());

		checkFolder(output_root / "imgs");

#pragma omp parallel for
		for (int i_img = 0; i_img < img_paths.size(); ++i_img)
		{
			const auto in_img_path = (fs::path(img_path) / img_paths[i_img]).string();
			const auto out_img_path = (output_root / "imgs" / img_paths[i_img]).string();
			cv::Mat img = cv::imread(in_img_path, cv::IMREAD_ANYCOLOR);
			cv::imwrite(out_img_path, img);

			Exiv2::Image::AutoPtr image1 = Exiv2::ImageFactory::open(in_img_path);
			assert(image1.get() != 0);
			image1->readMetadata();
			Exiv2::ExifData& exifData = image1->exifData();
			if (exifData.empty())
			{
				LOG(ERROR) << "No Exif data found in the file " << in_img_path;
				throw;
			}
			exifData["Exif.Image.Orientation"].setValue("1");
			exifData["Exif.Thumbnail.Orientation"].setValue("1");

			Exiv2::Image::AutoPtr image2 = Exiv2::ImageFactory::open(out_img_path);
			image2->setExifData(exifData);
			image2->writeMetadata();
		}
	}

	std::vector<std::string> ref_image_names;
	std::vector<Eigen::Vector3d> ref_locations;

	// Get gps location and proj
	LOG(INFO) << "===== Get gps location and proj =====";
	{
		fs::directory_iterator it_dir(output_root / "imgs"), it_dir_end;
		auto C = proj_context_create();

		// transform config
		const char* espg_4326 = "+proj=longlat +datum=WGS84 +no_defs ";
		const char* espg_3857 =
			"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 "
			"+y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs";
		const char* espg_4547 =
			"+proj=tmerc +lat_0=0 +lon_0=114 +k=1 +x_0=500000 +y_0=0 "
			"+ellps=GRS80 +units=m +no_defs";

		auto coordinate_transformer = proj_create_crs_to_crs(
			C,
			espg_4326,
			use_4547 ? espg_4547 : espg_3857,
			NULL);

		if (0 == coordinate_transformer) {
			fprintf(stderr, "Failed to create transformation object.\n");
			throw;
		}

		for (; it_dir != it_dir_end; ++it_dir)
		{
			ref_image_names.emplace_back(it_dir->path().filename().string());
			auto pos = getImageLocation_revised(it_dir->path().string(), coordinate_transformer, origin_point);
			ref_locations.emplace_back(pos);
		}

		proj_destroy(coordinate_transformer);
		proj_context_destroy(C); /* may be omitted in the single threaded case */
		assert(ref_image_names.size() == ref_locations.size());
	}

	//return 0;

	//Colmap running
	LOG(INFO) << "===== Colmap Running =====";
	ReconstructionManager reconstruction_manager;
	bool passAlignment = false;
	// Running Dense recon
	if (fs::exists(output_root / "sparse_align0")) {
		LOG(INFO) << "===== Existing Aligned Sparse Reconstruction, Dense Reconstruction Running =====";
		passAlignment = true;
		std::int16_t i = -1;
		// maybe sparse0,1,2.. exists
		while (fs::exists(output_root / ("sparse_align" + std::to_string(++i)))) {
			reconstruction_manager.Read((output_root / ("sparse_align" + std::to_string(i))).string());
		}
	}
	// Just running the alignment on [sparse0,1,2..] if exits
	else if (fs::exists(output_root / "sparse0")) {
		LOG(INFO) << "===== Existing Sparse Reconstruction, Alignment Only =====";
		std::int16_t i = -1;
		// maybe sparse0,1,2.. exists
		while (fs::exists(output_root / ("sparse" + std::to_string(++i)))) {
			reconstruction_manager.Read((output_root / ("sparse" + std::to_string(i))).string());
		}
	}
	// Recon
	else {
		LOG(INFO) << "===== Beginning Sparse Reconstruction =====";
		Database database(db_path);

		/*************************************** SIFT ***************************************/
		{
			ImageReaderOptions image_reader_options;
			SiftExtractionOptions sift_extraction_options;

			image_reader_options.database_path = db_path;
			image_reader_options.image_path = img_path;
			image_reader_options.camera_model = "SIMPLE_RADIAL";
			image_reader_options.single_camera = isSingleCamera;
			//image_reader_options.single_camera_per_folder   = true;

			sift_extraction_options.use_gpu = isUseGpu;
			//sift_extraction_options.upright = true;
			sift_extraction_options.gpu_index = "0";
			sift_extraction_options.max_image_size = max_image_size;

			if (isHighQuality)
			{
				//sift_extraction_options.estimate_affine_shape = true;
				sift_extraction_options.domain_size_pooling = true;
				sift_extraction_options.max_num_features = 40000;
				//sift_extraction_options.num_threads = 6;
			}

			SiftFeatureExtractor sift(image_reader_options, sift_extraction_options);
			sift.Start();
			sift.Wait();
			sift.Stop();
		}

		// set the pior
		{
			//set cameras prior(focal length)
			std::vector<Camera> cameras = database.ReadAllCameras();
			for (auto& camera : cameras) {
				camera.SetFocalLength(focal_length_in_pixel);
				camera.SetPriorFocalLength(true);
				database.UpdateCamera(camera);
			}
			// set image prior
			std::vector<Image> images = database.ReadAllImages();
			for (auto& image : images) {
				auto find_result = std::find(ref_image_names.begin(), ref_image_names.end(), image.Name());
				assert(find_result != ref_image_names.end());
				auto pos = ref_locations[find_result - ref_image_names.begin()];
				image.SetTvecPrior(pos);
				database.UpdateImage(image);
			}
		}

		// feature matching
		{
			SiftMatchingOptions sift_matching_options;
			ExhaustiveMatchingOptions exhaustive_matching_options;
			sift_matching_options.use_gpu = isUseGpu;
			sift_matching_options.gpu_index = "0";
			if (isHighQuality)
			{
				sift_matching_options.guided_matching = true;
				//sift_matching_options.planar_scene    = true;
				//sift_matching_options.max_num_matches = 43000;
			}
			ExhaustiveFeatureMatcher matcher(exhaustive_matching_options, sift_matching_options, db_path);
			matcher.Start();
			matcher.Wait();
			matcher.Stop();
		}

		/*************************************** SFM ***************************************/
		IncrementalMapperOptions incremental_mapper_options;

		// to fix intr and extr
		// // If reconstruction is provided as input, fix the existing image poses.
		// incremental_mapper_options.fix_existing_images = true;
		// // Which intrinsic parameters to optimize during the reconstruction.
		// incremental_mapper_options.ba_refine_focal_length = false;
		// incremental_mapper_options.ba_refine_extra_params = false;
		// incremental_mapper_options.ba_refine_principal_point = true;

		// recon option
		incremental_mapper_options.multiple_models = isMultipleModels;
		incremental_mapper_options.min_num_matches = 100; // 15 for more points +100

		// ba option
		if (isGreaterBa)
		{
			incremental_mapper_options.ba_local_max_num_iterations = 40;
			incremental_mapper_options.ba_local_max_refinements = 3;
			incremental_mapper_options.ba_global_max_num_iterations = 100;
			incremental_mapper_options.ba_global_use_pba = true;
		}

		IncrementalMapperController mapper(&incremental_mapper_options, img_path, db_path, &reconstruction_manager);
		mapper.Start();
		mapper.Wait();
		mapper.Stop();

		LOG(INFO) << boost::format("Reconstruction_manager.Size():%d") % (int)reconstruction_manager.Size();
		// for each recon model
		for (int i = 0; i < reconstruction_manager.Size(); i++)
		{
			LOG(INFO) << boost::format("========== Sparse Reconstruction %d/%d is Processing ==========") % (i + 1)
				% (int)reconstruction_manager.Size();

			Reconstruction& reconstruction = reconstruction_manager.Get(i);

			// Print some RegImageIds info
			for (const image_t i_img : reconstruction.RegImageIds()) {
				if (!reconstruction.IsImageRegistered(i_img))
					continue;
				LOG(INFO) << boost::format("Img %d: %d; ; Matched points: %d; Valid points: %d")
					% i_img
					% database.NumKeypointsForImage(i_img)
					% reconstruction.Image(i_img).NumCorrespondences()
					% reconstruction.Image(i_img).NumPoints3D();
			}
			LOG(INFO) << boost::format("%d registred images; %d 3D points")
				% reconstruction.NumRegImages()
				% reconstruction.NumPoints3D();

			// Export src result
			std::string outfile("sparse" + std::to_string(i));
			safeCheckFolder(output_root / outfile);
			reconstruction.WriteBinary((output_root / outfile).string());
			reconstruction.WriteText((output_root / outfile).string());
			//reconstruction.ExportPLY((output_root / outfile / "sparse.ply").string());
			ExportPLYwithNormal(reconstruction, (output_root / outfile / "sparse.ply").string());
			LOG(INFO) << boost::format("Reconstruction %d : Save to %s") % (i + 1) % (output_root / outfile);
		}
	}

	//return 0;
	// Alignment
	LOG(INFO) << boost::format("============= %d Reconstruction Begin Align =============\n") % (int)reconstruction_manager.Size();
	for (int i = 0; i < reconstruction_manager.Size(); i++)
	{
		// when folder ./sparse_align0 exists, it will be set to true
		if (passAlignment) {
			LOG(INFO) << boost::format("============= Passing Alignment =============\n");
			break;
		}

		LOG(INFO) << boost::format("============= Reconstruction %d/%d Begin Align =============\n") % (i + 1)
			% (int)reconstruction_manager.Size();

		Reconstruction& reconstruction = reconstruction_manager.Get(i);

		colmap::RANSACOptions ransac_options;
		SimilarityTransform3 tform;
		bool merge_origins = false;
		bool robust_alignment = true;
		int min_common_images = 3;
		bool alignment_success = true;

		// To find a best ransac_max_error
		std::double_t best_ransac_max_error = 0;
		std::vector<double> range = { 0.001,0.005,1.0 };
		std::vector<double> best_errors = { MAXINT,MAXINT,MAXINT,MAXINT }; // store the align error
		for (double ransac_max_error = range[0]; ransac_max_error <= range[2]; ransac_max_error += range[1])
		{
			std::cout << StringPrintf("=====> Testing ransac_max_error = %f <=====", ransac_max_error) << std::endl;
			ransac_options.max_error = ransac_max_error;
			//ransac_options.min_inlier_ratio             = 0.1;
			//ransac_options.dyn_num_trials_multiplier    = 100;

			alignment_success = reconstruction.AlignRobust(ref_image_names, ref_locations,
				min_common_images, ransac_options, &tform);

			if (!alignment_success)
			{
				std::cout << "Alignment failed\n" << std::endl;
				continue;
			}

			// alignment_success: error calculate
			std::vector<double> errors;
			errors.reserve(ref_image_names.size());
			for (size_t i = 0; i < ref_image_names.size(); ++i) {
				const Image* image = reconstruction.FindImageWithName(ref_image_names[i]);
				if (image != nullptr) {
					errors.push_back((image->ProjectionCenter() - ref_locations[i]).norm());
				}
				else {
					//LOG(INFO) << boost::format("Alignment: Image %s not used") % (ref_image_names[i]);
				}
			}

			std::cout << StringPrintf("Alignment error: %f (mean), %f (median), %f (max), %f (min)\n",
				Mean(errors), Median(errors), *std::max_element(errors.begin(), errors.end()), *std::min_element(errors.begin(), errors.end())) << std::endl;

			// if it's better, update the result
			std::double_t sumError = Mean(errors) + *std::max_element(errors.begin(), errors.end());
			if (Mean(errors) < best_errors[0]) {
				//std::cout << "Alignment: find better\n" << std::endl;
				best_errors[0] = Mean(errors);
				best_errors[1] = Median(errors);
				best_errors[2] = *std::max_element(errors.begin(), errors.end());
				best_errors[3] = *std::min_element(errors.begin(), errors.end());
				best_ransac_max_error = ransac_max_error;
			}
		}

		ransac_options.max_error = best_ransac_max_error;
		reconstruction.AlignRobust(ref_image_names, ref_locations, min_common_images, ransac_options, &tform);

		std::cout << "=====> Alignment End <=====" << std::endl;
		std::cout << StringPrintf("=====> Besr ransac_max_error = %f  <=====", best_ransac_max_error) << std::endl;
		std::cout << StringPrintf("Alignment error: %f (mean), %f (median), %f (max), %f (min)\n",
			best_errors[0], best_errors[1], best_errors[2], best_errors[3]) << std::endl;

		std::string outfile_align("sparse_align" + std::to_string(i));
		safeCheckFolder(output_root / outfile_align);
		reconstruction.WriteBinary((output_root / outfile_align).string());
		reconstruction.WriteText((output_root / outfile_align).string());
		reconstruction.ExportPLY((output_root / outfile_align / "sparse_align.ply").string());
		LOG(INFO) << boost::format("Aligned Reconstruction %d : Save to %s\n") % i % (output_root / outfile_align);

		/*int point3D_idx = 0;
		Point3D point3D = reconstruction.Point3D(point3D_idx);
		point3D.Error();*/
	}

	// Dense Reconstruction
	LOG(INFO) << "===== Beginning Dense Reconstruction =====";
	using namespace mvs;
	assert(reconstruction_manager.Size() != 0);
	for (int i = 0; i < reconstruction_manager.Size(); i++)
	{
		LOG(INFO) << boost::format("============= Dense Reconstruction %d/%d Begin =============\n") % (i + 1)
			% (int)reconstruction_manager.Size();
		std::string outfile_dense("dense" + std::to_string(i));
		std::string dense_path = (output_root / outfile_dense).string();
		std::string fused_path = (output_root / outfile_dense / "fused.ply").string();
		safeCheckFolder(output_root / outfile_dense);

		Reconstruction& reconstruction = reconstruction_manager.Get(i);

		// Undistort images and export undistorted cameras, as required by the mvs::PatchMatchController class.
		LOG(INFO) << "===== Beginning Images Undistortion =====";
		UndistortCameraOptions undistortion_options;
		undistortion_options.max_image_size = max_image_size;
		COLMAPUndistorter undistorter(undistortion_options, reconstruction, img_path, dense_path);
		undistorter.Start();
		undistorter.Wait();
		undistorter.Stop();

		// MVS
		LOG(INFO) << "===== Beginning Patch Match =====";
		PatchMatchOptions patchMatchOptions;
		patchMatchOptions.max_image_size = max_image_size;
		patchMatchOptions.gpu_index = "0";

		PatchMatchController patch_match_controller(patchMatchOptions, dense_path, "COLMAP", "");
		patch_match_controller.Start();
		patch_match_controller.Wait();
		patch_match_controller.Stop();

		// Stereo fusion
		LOG(INFO) << "===== Beginning Stereo Fusion =====";
		StereoFusionOptions stereoFusionOptions;
		stereoFusionOptions.max_image_size = max_image_size;
		const int num_reg_images = reconstruction.NumRegImages();
		stereoFusionOptions.min_num_pixels = std::min(num_reg_images + 1, stereoFusionOptions.min_num_pixels);
		StereoFusion fuser(stereoFusionOptions, dense_path, "COLMAP", "", "geometric");
		fuser.Start();
		fuser.Wait();
		fuser.Stop();

		std::cout << "Writing output: " << fused_path << std::endl;
		WriteBinaryPlyPoints(fused_path, fuser.GetFusedPoints());
		WritePointsVisibility(fused_path + ".vis", fuser.GetFusedPointsVisibility());
		LOG(INFO) << boost::format("Dense Reconstruction %d : Save to %s") % (i + 1) % (fused_path);

		// Surface meshing
		LOG(INFO) << "===== Beginning Surface Meshing =====";
		std::string meshing_path = (output_root / outfile_dense / "meshed-poisson.ply").string();
		PoissonMeshingOptions poissonMeshingOptions;
		poissonMeshingOptions.num_threads = 4;
		poissonMeshingOptions.trim = 1;
		PoissonMeshing(poissonMeshingOptions, fused_path, meshing_path);
		LOG(INFO) << boost::format("PoissonMesh %d : Save to %s") % (i + 1) % (meshing_path);
	}
	system("pause");
}

int runcolmap()
{
	fs::path output_root = "H:/data/runcolmapData";
	const std::string db_path = "H:/data/runcolmapData/database.db";  // Canbe blank
	const std::string img_path = "H:/data/runcolmapData/images";      // Path of the input image

	// Some config option
	// const double focal_length_in_pixel = 9004;

	// Recon option
	bool isHighQuality = false;
	bool isMultipleModels = false;
	bool isUseGpu = true;
	bool isSingleCamera = false;
	bool isGreaterBa = true;
	std::int16_t max_image_size = 6000;

	std::vector<std::string> ref_image_names;
	std::vector<Eigen::Vector3d> ref_locations;

	//Colmap running
	LOG(INFO) << "===== Colmap Running =====";
	ReconstructionManager reconstruction_manager;
	// sparse recon exists, loading
	if (fs::exists(output_root / "sparse0")) {
		LOG(INFO) << "===== Existing Aligned Sparse Reconstruction, Dense Reconstruction Running =====";
		std::int16_t i = -1;
		// maybe sparse0,1,2.. exists
		while (fs::exists(output_root / ("sparse" + std::to_string(++i)))) {
			reconstruction_manager.Read((output_root / ("sparse" + std::to_string(i))).string());
		}
	}
	// sparse recon
	else {
		LOG(INFO) << "===== Beginning Sparse Reconstruction =====";
		Database database(db_path);

		/*************************************** SIFT ***************************************/
		// SIFT config
		{
			ImageReaderOptions image_reader_options;
			SiftExtractionOptions sift_extraction_options;

			image_reader_options.database_path = db_path;
			image_reader_options.image_path = img_path;
			image_reader_options.camera_model = "SIMPLE_RADIAL";
			image_reader_options.single_camera = isSingleCamera;
			//image_reader_options.single_camera_per_folder   = true;

			sift_extraction_options.use_gpu = isUseGpu;
			//sift_extraction_options.upright = true;
			sift_extraction_options.gpu_index = "0";
			sift_extraction_options.max_image_size = max_image_size;

			if (isHighQuality)
			{
				//sift_extraction_options.estimate_affine_shape = true;
				sift_extraction_options.domain_size_pooling = true;
				sift_extraction_options.max_num_features = 40000;
				//sift_extraction_options.num_threads = 6;
			}

			SiftFeatureExtractor sift(image_reader_options, sift_extraction_options);
			sift.Start();
			sift.Wait();
			sift.Stop();
		}

		// set the pior
		/* {
			//set cameras prior(focal length)
			std::vector<Camera> cameras = database.ReadAllCameras();
			for (auto& camera : cameras) {
				camera.SetFocalLength(focal_length_in_pixel);
				camera.SetPriorFocalLength(true);
				database.UpdateCamera(camera);
			}
			// set image prior
			std::vector<Image> images = database.ReadAllImages();
			for (auto& image : images) {
				auto find_result = std::find(ref_image_names.begin(), ref_image_names.end(), image.Name());
				assert(find_result != ref_image_names.end());
				auto pos = ref_locations[find_result - ref_image_names.begin()];
				image.SetTvecPrior(pos);
				database.UpdateImage(image);
			}
		}*/

		// feature matching
		{
			SiftMatchingOptions sift_matching_options;
			ExhaustiveMatchingOptions exhaustive_matching_options;
			sift_matching_options.use_gpu = isUseGpu;
			sift_matching_options.gpu_index = "0";

			if (isHighQuality)
			{
				sift_matching_options.guided_matching = true;
				//sift_matching_options.planar_scene    = true;
				//sift_matching_options.max_num_matches = 43000;
			}
			ExhaustiveFeatureMatcher matcher(exhaustive_matching_options, sift_matching_options, db_path);
			matcher.Start();
			matcher.Wait();
			matcher.Stop();
		}

		/*************************************** SFM ***************************************/
		IncrementalMapperOptions incremental_mapper_options;

		// to fix intr and extr
		// // If reconstruction is provided as input, fix the existing image poses.
		// incremental_mapper_options.fix_existing_images = true;
		// // Which intrinsic parameters to optimize during the reconstruction.
		// incremental_mapper_options.ba_refine_focal_length = false;
		// incremental_mapper_options.ba_refine_extra_params = false;
		// incremental_mapper_options.ba_refine_principal_point = true;

		// recon option
		incremental_mapper_options.multiple_models = isMultipleModels;
		incremental_mapper_options.min_num_matches = 100; // 15 for more points +100

		// ba option
		if (isGreaterBa)
		{
			incremental_mapper_options.ba_local_max_num_iterations = 40;
			incremental_mapper_options.ba_local_max_refinements = 3;
			incremental_mapper_options.ba_global_max_num_iterations = 100;
			incremental_mapper_options.ba_global_use_pba = true;
		}

		IncrementalMapperController mapper(&incremental_mapper_options, img_path, db_path, &reconstruction_manager);
		mapper.Start();
		mapper.Wait();
		mapper.Stop();

		LOG(INFO) << boost::format("Reconstruction_manager.Size():%d") % (int)reconstruction_manager.Size();
		// for each recon model
		for (int i = 0; i < reconstruction_manager.Size(); i++)
		{
			LOG(INFO) << boost::format("========== Sparse Reconstruction %d/%d is Processing ==========") % (i + 1)
				% (int)reconstruction_manager.Size();

			Reconstruction& reconstruction = reconstruction_manager.Get(i);

			// Print some RegImageIds info
			for (const image_t i_img : reconstruction.RegImageIds()) {
				if (!reconstruction.IsImageRegistered(i_img))
					continue;
				LOG(INFO) << boost::format("Img %d: %d; ; Matched points: %d; Valid points: %d")
					% i_img
					% database.NumKeypointsForImage(i_img)
					% reconstruction.Image(i_img).NumCorrespondences()
					% reconstruction.Image(i_img).NumPoints3D();
			}
			LOG(INFO) << boost::format("%d registred images; %d 3D points")
				% reconstruction.NumRegImages()
				% reconstruction.NumPoints3D();

			// Export src result
			std::string outfile("sparse" + std::to_string(i));
			safeCheckFolder(output_root / outfile);
			reconstruction.WriteBinary((output_root / outfile).string());
			reconstruction.WriteText((output_root / outfile).string());
			//reconstruction.ExportPLY((output_root / outfile / "sparse.ply").string());
			ExportPLYwithNormal(reconstruction, (output_root / outfile / "sparse.ply").string());
			LOG(INFO) << boost::format("Reconstruction %d : Save to %s") % (i + 1) % (output_root / outfile);
		}
	}

	// dense Reconstruction
	LOG(INFO) << "===== Beginning Dense Reconstruction =====";
	using namespace mvs;
	assert(reconstruction_manager.Size() != 0);
	for (int i = 0; i < reconstruction_manager.Size(); i++)
	{
		LOG(INFO) << boost::format("============= Dense Reconstruction %d/%d Begin =============\n") % (i + 1)
			% (int)reconstruction_manager.Size();
		std::string outfile_dense("dense" + std::to_string(i));
		std::string dense_path = (output_root / outfile_dense).string();
		std::string fused_path = (output_root / outfile_dense / "fused.ply").string();
		safeCheckFolder(output_root / outfile_dense);

		Reconstruction& reconstruction = reconstruction_manager.Get(i);

		// Undistort images and export undistorted cameras, as required by the mvs::PatchMatchController class.
		LOG(INFO) << "===== Beginning Images Undistortion =====";
		UndistortCameraOptions undistortion_options;
		undistortion_options.max_image_size = max_image_size;
		COLMAPUndistorter undistorter(undistortion_options, reconstruction, img_path, dense_path);
		undistorter.Start();
		undistorter.Wait();
		undistorter.Stop();

		// MVS
		LOG(INFO) << "===== Beginning Patch Match =====";
		PatchMatchOptions patchMatchOptions;
		patchMatchOptions.max_image_size = max_image_size;
		patchMatchOptions.gpu_index = "0";
		patchMatchOptions.filter_min_triangulation_angle = 1.0f;

		PatchMatchController patch_match_controller(patchMatchOptions, dense_path, "COLMAP", "");
		patch_match_controller.Start();
		patch_match_controller.Wait();
		patch_match_controller.Stop();

		// Stereo fusion
		LOG(INFO) << "===== Beginning Stereo Fusion =====";
		StereoFusionOptions stereoFusionOptions;
		stereoFusionOptions.max_image_size = max_image_size;
		const int num_reg_images = reconstruction.NumRegImages();
		stereoFusionOptions.min_num_pixels = std::min(num_reg_images + 1, stereoFusionOptions.min_num_pixels);
		StereoFusion fuser(stereoFusionOptions, dense_path, "COLMAP", "", "geometric");
		fuser.Start();
		fuser.Wait();
		fuser.Stop();

		std::cout << "Writing output: " << fused_path << std::endl;
		WriteBinaryPlyPoints(fused_path, fuser.GetFusedPoints());
		WritePointsVisibility(fused_path + ".vis", fuser.GetFusedPointsVisibility());
		LOG(INFO) << boost::format("Dense Reconstruction %d : Save to %s") % (i + 1) % (fused_path);

		// Surface meshing
		LOG(INFO) << "===== Beginning Surface Meshing =====";
		std::string meshing_path = (output_root / outfile_dense / "meshed-poisson.ply").string();
		PoissonMeshingOptions poissonMeshingOptions;
		poissonMeshingOptions.num_threads = 4;
		poissonMeshingOptions.trim = 1;
		PoissonMeshing(poissonMeshingOptions, fused_path, meshing_path);
		LOG(INFO) << boost::format("PoissonMesh %d : Save to %s") % (i + 1) % (meshing_path);
	}
}

// std::string scene_src_path, std::string scene_output_path
int recon_scene()
{
	const fs::path input_root("H:/data/temp/scene0000_00");
	const fs::path output_root("H:/data/scannet_colmap/scene0000_00");

	const std::string img_path = (input_root / "color").string();			// Path of the input image
	const fs::path pose_folder_path = input_root / "pose";					// c2w pose
	const fs::path intr_txt_path = input_root / "intrinsics_color.txt";		// c2w pose
	const std::string db_path = (input_root / "database.db").string();		// Canbe blank

	// Some config option
	// Recon option
	bool isHighQuality = false;
	bool isMultipleModels = false;
	bool isUseGpu = true;
	bool isSingleCamera = false;
	bool isGreaterBa = true;
	std::int16_t max_image_size = 6000;

	//Colmap running
	LOG(INFO) << "===== Colmap Running =====";
	ReconstructionManager reconstruction_manager;
	// sparse recon exists, loading
	if (fs::exists(output_root / "sparse0")) {
		LOG(INFO) << "===== Existing Aligned Sparse Reconstruction, Dense Reconstruction Running =====";
		std::int16_t i = -1;
		// maybe sparse0,1,2.. exists
		while (fs::exists(output_root / ("sparse" + std::to_string(++i)))) {
			reconstruction_manager.Read((output_root / ("sparse" + std::to_string(i))).string());
		}
	}
	// sparse recon
	else {
		LOG(INFO) << "===== Beginning Sparse Reconstruction =====";
		Database database(db_path);

		/*************************************** SIFT ***************************************/
		// SIFT config
		{
			ImageReaderOptions image_reader_options;
			SiftExtractionOptions sift_extraction_options;

			image_reader_options.database_path = db_path;
			image_reader_options.image_path = img_path;
			image_reader_options.camera_model = "SIMPLE_RADIAL";
			image_reader_options.single_camera = isSingleCamera;
			//image_reader_options.single_camera_per_folder   = true;

			sift_extraction_options.use_gpu = isUseGpu;
			//sift_extraction_options.upright = true;
			sift_extraction_options.gpu_index = "0";
			sift_extraction_options.max_image_size = max_image_size;

			if (isHighQuality)
			{
				//sift_extraction_options.estimate_affine_shape = true;
				sift_extraction_options.domain_size_pooling = true;
				sift_extraction_options.max_num_features = 40000;
				//sift_extraction_options.num_threads = 6;
			}

			SiftFeatureExtractor sift(image_reader_options, sift_extraction_options);
			sift.Start();
			sift.Wait();
			sift.Stop();
		}

		// feature matching
		{
			SiftMatchingOptions sift_matching_options;
			ExhaustiveMatchingOptions exhaustive_matching_options;
			sift_matching_options.use_gpu = isUseGpu;
			sift_matching_options.gpu_index = "0";

			if (isHighQuality)
			{
				sift_matching_options.guided_matching = true;
				//sift_matching_options.planar_scene    = true;
				//sift_matching_options.max_num_matches = 43000;
			}
			ExhaustiveFeatureMatcher matcher(exhaustive_matching_options, sift_matching_options, db_path);
			matcher.Start();
			matcher.Wait();
			matcher.Stop();
		}

		/********************************** update database *********************************/
		// set the pior
		{
			// parse the intr_txt and get the focal length
			Eigen::Matrix<double, 4, 4> intr_matrix = ReadMatrixFromFile(intr_txt_path.string());
			const double focal_length_in_pixel = intr_matrix(0, 0);

			//set cameras prior(focal length)
			std::vector<Camera> cameras = database.ReadAllCameras();
			for (auto& camera : cameras) {
				camera.SetFocalLength(focal_length_in_pixel);
				camera.SetPriorFocalLength(true);
				database.UpdateCamera(camera);
			}

			// for each image in database, we set and fixed its pose
			std::vector<Image> images = database.ReadAllImages();
			for (auto& image : images) {
				std::string txtName = image.Name();
				size_t pos = txtName.find_last_of(".");
				if (pos != std::string::npos) {
					txtName.replace(pos + 1, txtName.length() - pos - 1, "txt");
				}

				Eigen::Matrix<double, 4, 4> extr_matrix = ReadMatrixFromFile((pose_folder_path / txtName).string()).inverse();
				Eigen::Quaterniond q(extr_matrix.block<3, 3>(0, 0));

				image.SetQvec(Eigen::Vector4d(q.w(), q.x(), q.y(), q.z()));
				image.SetTvec(extr_matrix.col(3).head<3>());

				image.SetQvecPrior(Eigen::Vector4d(q.w(), q.x(), q.y(), q.z()));
				image.SetTvecPrior(extr_matrix.col(3).head<3>());
				database.UpdateImage(image);
			}
		}

		/*************************************** SFM ***************************************/
		IncrementalMapperOptions incremental_mapper_options;

		// fix the image pose
		incremental_mapper_options.fix_existing_images = true;
		// Which intrinsic parameters to optimize during the reconstruction.
		incremental_mapper_options.ba_refine_focal_length = false;
		incremental_mapper_options.ba_refine_extra_params = false;
		incremental_mapper_options.ba_refine_principal_point = true;

		// recon option
		incremental_mapper_options.multiple_models = isMultipleModels;
		incremental_mapper_options.min_num_matches = 100; // 15 for more points +100

		// ba option
		if (isGreaterBa)
		{
			incremental_mapper_options.ba_local_max_num_iterations = 40;
			incremental_mapper_options.ba_local_max_refinements = 3;
			incremental_mapper_options.ba_global_max_num_iterations = 100;
			incremental_mapper_options.ba_global_use_pba = true;
		}

		IncrementalMapperController mapper(&incremental_mapper_options, img_path, db_path, &reconstruction_manager);
		mapper.Start();
		mapper.Wait();
		mapper.Stop();

		LOG(INFO) << boost::format("Reconstruction_manager.Size():%d") % (int)reconstruction_manager.Size();
		// for each recon model
		for (int i = 0; i < reconstruction_manager.Size(); i++)
		{
			LOG(INFO) << boost::format("========== Sparse Reconstruction %d/%d is Processing ==========") % (i + 1)
				% (int)reconstruction_manager.Size();

			Reconstruction& reconstruction = reconstruction_manager.Get(i);

			// Print some RegImageIds info
			for (const image_t i_img : reconstruction.RegImageIds()) {
				if (!reconstruction.IsImageRegistered(i_img))
					continue;
				LOG(INFO) << boost::format("Img %d: %d; ; Matched points: %d; Valid points: %d")
					% i_img
					% database.NumKeypointsForImage(i_img)
					% reconstruction.Image(i_img).NumCorrespondences()
					% reconstruction.Image(i_img).NumPoints3D();
			}
			LOG(INFO) << boost::format("%d registred images; %d 3D points")
				% reconstruction.NumRegImages()
				% reconstruction.NumPoints3D();

			// Export src result
			std::string outfile("sparse" + std::to_string(i));
			safeCheckFolder(output_root / outfile);
			reconstruction.WriteBinary((output_root / outfile).string());
			reconstruction.WriteText((output_root / outfile).string());
			//reconstruction.ExportPLY((output_root / outfile / "sparse.ply").string());
			ExportPLYwithNormal(reconstruction, (output_root / outfile / "sparse.ply").string());
			LOG(INFO) << boost::format("Reconstruction %d : Save to %s") % (i + 1) % (output_root / outfile);
		}
	}

	// dense Reconstruction
	LOG(INFO) << "===== Beginning Dense Reconstruction =====";
	using namespace mvs;
	assert(reconstruction_manager.Size() != 0);
	for (int i = 0; i < reconstruction_manager.Size(); i++)
	{
		LOG(INFO) << boost::format("============= Dense Reconstruction %d/%d Begin =============\n") % (i + 1)
			% (int)reconstruction_manager.Size();
		std::string outfile_dense("dense" + std::to_string(i));
		std::string dense_path = (output_root / outfile_dense).string();
		std::string fused_path = (output_root / outfile_dense / "fused.ply").string();
		safeCheckFolder(output_root / outfile_dense);

		Reconstruction& reconstruction = reconstruction_manager.Get(i);

		// Undistort images and export undistorted cameras, as required by the mvs::PatchMatchController class.
		LOG(INFO) << "===== Beginning Images Undistortion =====";
		UndistortCameraOptions undistortion_options;
		undistortion_options.max_image_size = max_image_size;
		COLMAPUndistorter undistorter(undistortion_options, reconstruction, img_path, dense_path);
		undistorter.Start();
		undistorter.Wait();
		undistorter.Stop();

		// MVS
		LOG(INFO) << "===== Beginning Patch Match =====";
		PatchMatchOptions patchMatchOptions;
		patchMatchOptions.max_image_size = max_image_size;
		patchMatchOptions.gpu_index = "0";
		// patchMatchOptions.filter_min_triangulation_angle = 1.0f;

		PatchMatchController patch_match_controller(patchMatchOptions, dense_path, "COLMAP", "");
		patch_match_controller.Start();
		patch_match_controller.Wait();
		patch_match_controller.Stop();

		// Stereo fusion
		LOG(INFO) << "===== Beginning Stereo Fusion =====";
		StereoFusionOptions stereoFusionOptions;
		stereoFusionOptions.max_image_size = max_image_size;
		const int num_reg_images = reconstruction.NumRegImages();
		stereoFusionOptions.min_num_pixels = std::min(num_reg_images + 1, stereoFusionOptions.min_num_pixels);
		StereoFusion fuser(stereoFusionOptions, dense_path, "COLMAP", "", "geometric");
		fuser.Start();
		fuser.Wait();
		fuser.Stop();

		std::cout << "Writing output: " << fused_path << std::endl;
		WriteBinaryPlyPoints(fused_path, fuser.GetFusedPoints());
		WritePointsVisibility(fused_path + ".vis", fuser.GetFusedPointsVisibility());
		LOG(INFO) << boost::format("Dense Reconstruction %d : Save to %s") % (i + 1) % (fused_path);

		// Surface meshing
		LOG(INFO) << "===== Beginning Surface Meshing =====";
		std::string meshing_path = (output_root / outfile_dense / "meshed-poisson.ply").string();
		PoissonMeshingOptions poissonMeshingOptions;
		poissonMeshingOptions.num_threads = 4;
		poissonMeshingOptions.trim = 1;
		PoissonMeshing(poissonMeshingOptions, fused_path, meshing_path);
		LOG(INFO) << boost::format("PoissonMesh %d : Save to %s") % (i + 1) % (meshing_path);
	}
}

int main()
{
	runcolmap_with_align();
}