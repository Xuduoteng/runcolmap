#include "align.h"
#include <vector>
#include "colmap/exe/model.h"
#include "colmap/base/gps.h"
#include "colmap/base/pose.h"
#include "colmap/base/similarity_transform.h"
#include "colmap/estimators/coordinate_frame.h"
#include "colmap/util/misc.h"
#include "colmap/util/option_manager.h"
#include "colmap/util/threading.h"
#include "colmap/optim/ransac.h"

using namespace colmap;

std::vector<Eigen::Vector3d> ConvertCameraLocations(
    const bool ref_is_gps, const std::string& alignment_type,
    const std::vector<Eigen::Vector3d>& ref_locations) {
  if (ref_is_gps) {
    const GPSTransform gps_transform(GPSTransform::WGS84);
    if (alignment_type != "enu") {
      std::cout << "\nConverting Alignment Coordinates from GPS (lat/lon/alt) "
                   "to ECEF.\n";
      return gps_transform.EllToXYZ(ref_locations);
    } else {
      std::cout << "\nConverting Alignment Coordinates from GPS (lat/lon/alt) "
                   "to ENU.\n";
      return gps_transform.EllToENU(ref_locations, ref_locations[0](0),
                                    ref_locations[0](1));
    }
  } else {
    std::cout << "\nCartesian Alignment Coordinates extracted (MUST NOT BE "
                 "GPS coords!).\n";
    return ref_locations;
  }
}

void ReadFileCameraLocations(const std::string& ref_images_path,
                             const bool ref_is_gps,
                             const std::string& alignment_type,
                             std::vector<std::string>* ref_image_names,
                             std::vector<Eigen::Vector3d>* ref_locations) {
  for (const auto& line : ReadTextFileLines(ref_images_path)) {
    std::stringstream line_parser(line);
    std::string image_name;
    Eigen::Vector3d camera_position;
    line_parser >> image_name >> camera_position[0] >> camera_position[1] >>
        camera_position[2];
    ref_image_names->push_back(image_name);
    ref_locations->push_back(camera_position);
  }

  *ref_locations =
      ConvertCameraLocations(ref_is_gps, alignment_type, *ref_locations);
}

void ReadDatabaseCameraLocations(const std::string& database_path,
                                 const bool ref_is_gps,
                                 const std::string& alignment_type,
                                 std::vector<std::string>* ref_image_names,
                                 std::vector<Eigen::Vector3d>* ref_locations) {
  Database database(database_path);
  for (const auto& image : database.ReadAllImages()) {
    if (image.HasTvecPrior()) {
      ref_image_names->push_back(image.Name());
      ref_locations->push_back(image.TvecPrior());
    }
  }

  *ref_locations =
      ConvertCameraLocations(ref_is_gps, alignment_type, *ref_locations);
}

// namespace

// Align given reconstruction with user provided cameras positions
// (can be used for geo-registration for instance).
// The cameras positions to be used for aligning the reconstruction
// model must be provided either by a txt file (with each line being: img_name x
// y z) or through a colmap database file containing a prior position for the
// registered images.
//
// Required Options:
// - input_path: path to initial reconstruction model
// - output_path: path to store the aligned reconstruction model
//
// Additional Options:
// - database_path: path to database file with prior positions for
// reconstruction images
// - ref_images_path: path to txt file with prior positions for reconstruction
// images (WARNING: provide only one of the above)
// - ref_is_gps: if true the prior positions are converted from GPS
// (lat/lon/alt) to ECEF or ENU
// - merge_image_and_ref_origins: if true the reconstuction will be shifted so
// that the first prior position is used for its camera position
// - transform_path: path to store the Sim3 transformation used for the
// alignment
// - alignment_type:
//    > plane: align with reconstruction principal plane
//    > ecef: align with ecef coords. (requires gps coords. or user provided
//    ecef coords.)
//    > enu: align with enu coords. (requires gps coords. or user provided enu
//    coords.)
//    > enu-plane: align to ecef and then to enu plane (requires gps
//    coords. or user provided ecef coords.)
//    > enu-plane-unscaled: same as above but do not apply the computed
//    scale when aligning the reconstruction
//    > custom: align to provided coords.
// - min_common_images: minimum number of images with prior positions to perform
// the estimate an alignment
// - estimate_scale: if true apply the computed scale when aligning the
// reconstruction
// - robust_alignment: if true use a ransac-based estimation for robust
// alignment
// - robust_alignment_max_error: ransac error to use if robust alignment is
// enabled
int align(std::string input_path,
                    std::string output_path,
                    //std::string database_path,
                    //std::string ref_images_path,
                    std::vector<std::string> &ref_image_names,
                    std::vector<Eigen::Vector3d> &ref_locations,
                    bool ref_is_gps,
                    bool merge_origins,
                    std::string alignment_type,
                    int min_common_images,
                    bool robust_alignment,
                    bool estimate_scale,
                    std::double_t max_error) {

  colmap::RANSACOptions ransac_options;
  ransac_options.max_error = max_error;
  std::string transform_path;

  colmap::StringToLower(&alignment_type);
  const std::unordered_set<std::string> alignment_options{
      "plane", "ecef", "enu", "enu-plane", "enu-plane-unscaled", "custom"};
  if (alignment_options.count(alignment_type) == 0) {
    std::cerr << "ERROR: Invalid `alignment_type` - supported values are "
                 "{'plane', 'ecef', 'enu', 'enu-plane', 'enu-plane-unscaled', "
                 "'custom'}"
              << std::endl;
    return EXIT_FAILURE;
  }

  if (robust_alignment && ransac_options.max_error <= 0) {
    std::cout << "ERROR: You must provide a maximum alignment error > 0"
              << std::endl;
    return EXIT_FAILURE;
  }

  // if (alignment_type != "plane" && database_path.empty() &&
  //     ref_images_path.empty()) {
  //   std::cerr << "ERROR: Location alignment requires either database or "
  //                "location file path."
  //             << std::endl;
  //   return EXIT_FAILURE;
  // }

  // 2
  // std::vector<std::string> ref_image_names;
  // std::vector<Eigen::Vector3d> ref_locations;
  // if (!ref_images_path.empty() && database_path.empty()) {
  //   std::cout << "ReadFileCameraLocations running" << std::endl;
  //   ReadFileCameraLocations(ref_images_path, ref_is_gps, alignment_type,
  //                           &ref_image_names, &ref_locations);
  // } else if (!database_path.empty() && ref_images_path.empty()) {
  //   std::cout << "ReadDatabaseCameraLocations running" << std::endl;
  //   ReadDatabaseCameraLocations(database_path, ref_is_gps, alignment_type,
  //                               &ref_image_names, &ref_locations);
  // } else if (alignment_type != "plane") {
  //   std::cerr << "ERROR: Use location file or database, not both" << std::endl;
  //   return EXIT_FAILURE;
  // }

  if (alignment_type != "plane" &&
      static_cast<int>(ref_locations.size()) < min_common_images) {
    std::cout << "ERROR: Cannot align with insufficient reference locations:" 
              << static_cast<int>(ref_locations.size())
              << std::endl;
    return EXIT_FAILURE;
  }

  Reconstruction reconstruction;
  reconstruction.Read(input_path);
  SimilarityTransform3 tform;
  bool alignment_success = true;

  if (alignment_type == "plane") {
    PrintHeading2("Aligning reconstruction to principal plane");
    AlignToPrincipalPlane(&reconstruction, &tform);
  } else {
    PrintHeading2("Aligning reconstruction to " + alignment_type);
    std::cout << StringPrintf(" => Using %d reference images",
                              ref_image_names.size())
              << std::endl;

    if (estimate_scale) {
      if (robust_alignment) {
        alignment_success = reconstruction.AlignRobust(
            ref_image_names, ref_locations, min_common_images, ransac_options,
            &tform);
      } else {
        alignment_success = reconstruction.Align(ref_image_names, ref_locations,
                                                 min_common_images, &tform);
      }
    } else {
      if (robust_alignment) {
        alignment_success = reconstruction.AlignRobust<false>(
            ref_image_names, ref_locations, min_common_images, ransac_options,
            &tform);
      } else {
        alignment_success = reconstruction.Align<false>(
            ref_image_names, ref_locations, min_common_images, &tform);
      }
    }

    std::vector<double> errors;
    errors.reserve(ref_image_names.size());
  
    for (size_t i = 0; i < ref_image_names.size(); ++i) {
      const Image* image = reconstruction.FindImageWithName(ref_image_names[i]);
      if (image != nullptr) {
        errors.push_back((image->ProjectionCenter() - ref_locations[i]).norm());
      }
      else {
        std::cout << "image == nullptr" << std::endl;
      }
    }

    std::cout << StringPrintf("=> Alignment error: %f (mean), %f (median)",
                              Mean(errors), Median(errors))
              << std::endl;


    if (alignment_success && StringStartsWith(alignment_type, "enu-plane")) {
      PrintHeading2("Aligning ECEF aligned reconstruction to ENU plane");
      AlignToENUPlane(&reconstruction, &tform,
                      alignment_type == "enu-plane-unscaled");
    }
  }

  if (merge_origins) {
    for (size_t i = 0; i < ref_image_names.size(); i++) {
      const Image* first_image =
          reconstruction.FindImageWithName(ref_image_names[i]);

      if (first_image != nullptr) {
        const Eigen::Vector3d& first_img_position = ref_locations[i];

        const Eigen::Vector3d trans_align =
            first_img_position - first_image->ProjectionCenter();

        const SimilarityTransform3 origin_align(
            1.0, ComposeIdentityQuaternion(), trans_align);

        std::cout << "\n Aligning Reconstruction's origin with Ref origin : "
                  << first_img_position.transpose() << "\n";

        reconstruction.Transform(origin_align);

        // Update the Sim3 transformation in case it is stored next
        tform = SimilarityTransform3(tform.Scale(), tform.Rotation(),
                                     tform.Translation() + trans_align);

        break;
      }
    }
  }

  if (alignment_success) {
    std::cout << "=> Alignment succeeded" << std::endl;
    reconstruction.Write(output_path);
    reconstruction.WriteText(output_path);
    //reconstruction.ExportPLY(output_path + "/sparse_align.ply");
    if (!transform_path.empty()) {
      tform.Write(transform_path);
    }
    return EXIT_SUCCESS;
  } else {
    std::cout << "=> Alignment failed" << std::endl;
    return EXIT_FAILURE;
  }
}