#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

std::vector<Eigen::Vector3d> ConvertCameraLocations(
    const bool ref_is_gps, const std::string& alignment_type,
    const std::vector<Eigen::Vector3d>& ref_locations);

void ReadFileCameraLocations(const std::string& ref_images_path,
                             const bool ref_is_gps,
                             const std::string& alignment_type,
                             std::vector<std::string>* ref_image_names,
                             std::vector<Eigen::Vector3d>* ref_locations);


void ReadDatabaseCameraLocations(const std::string& database_path,
                                 const bool ref_is_gps,
                                 const std::string& alignment_type,
                                 std::vector<std::string>* ref_image_names,
                                 std::vector<Eigen::Vector3d>* ref_locations);


int align(std::string input_path,
                    std::string output_path,
                    //std::string database_path,
                    //std::string ref_images_path,
                    std::vector<std::string> &ref_image_names,
                    std::vector<Eigen::Vector3d> &ref_locations,
                    bool ref_is_gps = true,
                    bool merge_origins = false,
                    std::string alignment_type = "ecef",
                    int min_common_images = 3,
                    bool robust_alignment = true,
                    bool estimate_scale = true,
                    std::double_t max_error = 0.01);