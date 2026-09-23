#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace
{

constexpr double kPi = 3.141592653589793238462643383279502884;

struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b)
{
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(const Vec3& a, double value)
{
  return {a.x * value, a.y * value, a.z * value};
}

double dot(const Vec3& a, const Vec3& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b)
{
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

Vec3 normalized(const Vec3& value)
{
  const double norm = std::sqrt(dot(value, value));
  if(!std::isfinite(norm) || norm <= std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Cannot normalize a degenerate vector");
  return value * (1.0 / norm);
}

struct Options
{
  std::string output_path;
  std::uint64_t seed = 20260904;
  std::size_t views = 10;
  std::size_t image_width = 96;
  std::size_t image_height = 72;
  double vertical_fov_degrees = 50.0;
  double voxel_size = 0.018;
  double base_depth_noise_sigma = 0.0012;
  double depth_quantization = 0.0005;
  double dropout_probability = 0.06;
};

struct Sample
{
  Vec3 noisy;
  Vec3 ground_truth;
  std::size_t observations = 1;
};

bool parsePositiveDouble(const std::string& text, double& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size() && std::isfinite(value) && value > 0.0;
  }
  catch(const std::exception&)
  {
    return false;
  }
}

bool parseProbability(const std::string& text, double& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    if(parsed != text.size() || !std::isfinite(value))
      return false;
  }
  catch(const std::exception&)
  {
    return false;
  }
  return value >= 0.0 && value < 1.0;
}

bool parseSize(const std::string& text, std::size_t& value)
{
  try
  {
    std::size_t parsed = 0;
    value = static_cast<std::size_t>(std::stoull(text, &parsed));
    return parsed == text.size();
  }
  catch(const std::exception&)
  {
    return false;
  }
}

bool parseSeed(const std::string& text, std::uint64_t& value)
{
  try
  {
    std::size_t parsed = 0;
    value = static_cast<std::uint64_t>(std::stoull(text, &parsed));
    return parsed == text.size();
  }
  catch(const std::exception&)
  {
    return false;
  }
}

void printUsage(const char* program)
{
  std::cerr
      << "Usage: " << program << " OUTPUT.xyz [options]\n"
      << "All coordinates and distance parameters use meters (m).\n"
      << "  --seed INTEGER --views COUNT --width PIXELS --height PIXELS\n"
      << "  --vertical-fov-deg VALUE --voxel-size VALUE_METERS\n"
      << "  --noise-sigma VALUE_METERS --depth-quantization VALUE_METERS\n"
      << "  --dropout-probability VALUE_IN_0_1\n";
}

bool parseOptions(int argc, char** argv, Options& options)
{
  if(argc < 2)
    return false;
  options.output_path = argv[1];
  for(int i = 2; i < argc; ++i)
  {
    if(i + 1 >= argc)
      return false;
    const std::string option = argv[i];
    const std::string value = argv[++i];
    if(option == "--seed")
    {
      if(!parseSeed(value, options.seed))
        return false;
    }
    else if(option == "--views")
    {
      if(!parseSize(value, options.views) || options.views < 6)
        return false;
    }
    else if(option == "--width")
    {
      if(!parseSize(value, options.image_width) || options.image_width < 8)
        return false;
    }
    else if(option == "--height")
    {
      if(!parseSize(value, options.image_height) || options.image_height < 8)
        return false;
    }
    else if(option == "--vertical-fov-deg")
    {
      if(!parsePositiveDouble(value, options.vertical_fov_degrees) ||
         options.vertical_fov_degrees >= 120.0)
        return false;
    }
    else if(option == "--voxel-size")
    {
      if(!parsePositiveDouble(value, options.voxel_size))
        return false;
    }
    else if(option == "--noise-sigma")
    {
      if(!parsePositiveDouble(value, options.base_depth_noise_sigma))
        return false;
    }
    else if(option == "--depth-quantization")
    {
      if(!parsePositiveDouble(value, options.depth_quantization))
        return false;
    }
    else if(option == "--dropout-probability")
    {
      if(!parseProbability(value, options.dropout_probability))
        return false;
    }
    else
    {
      return false;
    }
  }
  return true;
}

bool rayEllipsoidIntersection(const Vec3& origin,
                              const Vec3& direction,
                              const Vec3& center,
                              const Vec3& radii,
                              double& distance,
                              Vec3& point,
                              Vec3& normal)
{
  const Vec3 relative = origin - center;
  const double rx2 = radii.x * radii.x;
  const double ry2 = radii.y * radii.y;
  const double rz2 = radii.z * radii.z;
  const double a = direction.x * direction.x / rx2 +
                   direction.y * direction.y / ry2 +
                   direction.z * direction.z / rz2;
  const double b = 2.0 * (relative.x * direction.x / rx2 +
                          relative.y * direction.y / ry2 +
                          relative.z * direction.z / rz2);
  const double c = relative.x * relative.x / rx2 +
                   relative.y * relative.y / ry2 +
                   relative.z * relative.z / rz2 - 1.0;
  const double discriminant = b * b - 4.0 * a * c;
  if(discriminant < 0.0)
    return false;
  const double root = std::sqrt(discriminant);
  const double near_distance = (-b - root) / (2.0 * a);
  const double far_distance = (-b + root) / (2.0 * a);
  distance = near_distance > 0.0 ? near_distance : far_distance;
  if(distance <= 0.0)
    return false;
  point = origin + direction * distance;
  normal = normalized({(point.x - center.x) / rx2,
                       (point.y - center.y) / ry2,
                       (point.z - center.z) / rz2});
  return true;
}

fs::path companionPath(const fs::path& output, const std::string& suffix)
{
  return output.parent_path() /
         (output.stem().string() + suffix + output.extension().string());
}

void writePoint(std::ostream& out, const Vec3& point)
{
  out << point.x << ' ' << point.y << ' ' << point.z << '\n';
}

}  // namespace

int main(int argc, char** argv)
{
  Options options;
  if(!parseOptions(argc, argv, options))
  {
    printUsage(argv[0]);
    return 1;
  }

  try
  {
    const Vec3 center{0.0, 0.0, 0.75};
    const Vec3 radii{0.18, 0.13, 0.24};
    const std::size_t ring_views = options.views - 2;
    const double aspect = static_cast<double>(options.image_width) /
                          static_cast<double>(options.image_height);
    const double tan_half_fov = std::tan(
        0.5 * options.vertical_fov_degrees * kPi / 180.0);

    std::mt19937_64 random(options.seed);
    std::uniform_real_distribution<double> unit_random(0.0, 1.0);
    std::normal_distribution<double> standard_normal(0.0, 1.0);
    std::map<std::array<int, 3>, Sample> voxels;
    std::size_t cast_rays = 0;
    std::size_t surface_hits = 0;
    std::size_t accepted_measurements = 0;

    for(std::size_t view = 0; view < options.views; ++view)
    {
      Vec3 camera;
      if(view < ring_views)
      {
        const double angle = 2.0 * kPi * static_cast<double>(view) /
                             static_cast<double>(ring_views);
        const double height_offset = view % 2 == 0 ? 0.15 : -0.12;
        camera = {0.85 * std::cos(angle), 0.85 * std::sin(angle),
                  center.z + height_offset};
      }
      else
      {
        camera = view == ring_views
                     ? Vec3{0.08, 0.0, center.z + 0.85}
                     : Vec3{-0.08, 0.0, center.z - 0.85};
      }

      const Vec3 forward = normalized(center - camera);
      const Vec3 up_hint = std::abs(forward.z) > 0.9
                               ? Vec3{0.0, 1.0, 0.0}
                               : Vec3{0.0, 0.0, 1.0};
      const Vec3 right = normalized(cross(forward, up_hint));
      const Vec3 up = normalized(cross(right, forward));

      for(std::size_t row = 0; row < options.image_height; ++row)
      {
        for(std::size_t column = 0; column < options.image_width; ++column)
        {
          ++cast_rays;
          const double image_x =
              (2.0 * (static_cast<double>(column) + 0.5) /
                   static_cast<double>(options.image_width) -
               1.0) *
              aspect * tan_half_fov;
          const double image_y =
              (1.0 - 2.0 * (static_cast<double>(row) + 0.5) /
                         static_cast<double>(options.image_height)) *
              tan_half_fov;
          const Vec3 direction = normalized(
              forward + right * image_x + up * image_y);

          double true_depth = 0.0;
          Vec3 true_point;
          Vec3 normal;
          if(!rayEllipsoidIntersection(camera, direction, center, radii,
                                       true_depth, true_point, normal))
            continue;
          ++surface_hits;

          const double incidence =
              std::clamp(-dot(direction, normal), 0.0, 1.0);
          const double dropout = std::min(
              0.8, options.dropout_probability + 0.22 * (1.0 - incidence));
          if(unit_random(random) < dropout)
            continue;

          const double depth_sigma = options.base_depth_noise_sigma *
                                     (1.0 + 1.5 * (1.0 - incidence));
          double noisy_depth =
              true_depth + depth_sigma * standard_normal(random);
          noisy_depth = options.depth_quantization *
                        std::round(noisy_depth / options.depth_quantization);
          const double tangential_sigma = 0.25 * depth_sigma;
          const Vec3 noisy_point =
              camera + direction * noisy_depth +
              right * (tangential_sigma * standard_normal(random)) +
              up * (tangential_sigma * standard_normal(random));
          ++accepted_measurements;

          const std::array<int, 3> key{
              static_cast<int>(std::floor(noisy_point.x /
                                          options.voxel_size)),
              static_cast<int>(std::floor(noisy_point.y /
                                          options.voxel_size)),
              static_cast<int>(std::floor(noisy_point.z /
                                          options.voxel_size))};
          auto inserted = voxels.emplace(
              key, Sample{noisy_point, true_point, 1});
          if(!inserted.second)
          {
            Sample& sample = inserted.first->second;
            ++sample.observations;
            // Reservoir sampling prevents the first camera in the loop from
            // systematically owning every multiply observed voxel.
            if(unit_random(random) <
               1.0 / static_cast<double>(sample.observations))
            {
              sample.noisy = noisy_point;
              sample.ground_truth = true_point;
            }
          }
        }
      }
    }

    if(voxels.empty())
      throw std::runtime_error("The simulated cameras produced no points");

    double noise_sum = 0.0;
    double noise_squared_sum = 0.0;
    double maximum_noise = 0.0;
    for(const auto& entry : voxels)
    {
      const Vec3 error = entry.second.noisy - entry.second.ground_truth;
      const double magnitude = std::sqrt(dot(error, error));
      noise_sum += magnitude;
      noise_squared_sum += magnitude * magnitude;
      maximum_noise = std::max(maximum_noise, magnitude);
    }
    const double mean_noise = noise_sum / static_cast<double>(voxels.size());
    const double rms_noise =
        std::sqrt(noise_squared_sum / static_cast<double>(voxels.size()));

    const fs::path output(options.output_path);
    if(!output.parent_path().empty())
      fs::create_directories(output.parent_path());
    const fs::path ground_truth_path =
        companionPath(output, "_ground_truth");
    const fs::path metadata_path =
        output.parent_path() / (output.stem().string() + "_metadata.json");
    std::ofstream noisy_output(output);
    std::ofstream ground_truth_output(ground_truth_path);
    if(!noisy_output || !ground_truth_output)
      throw std::runtime_error("Failed to open point-cloud output files");
    noisy_output << std::setprecision(17);
    ground_truth_output << std::setprecision(17);
    for(const auto& entry : voxels)
    {
      writePoint(noisy_output, entry.second.noisy);
      writePoint(ground_truth_output, entry.second.ground_truth);
    }

    std::ofstream metadata(metadata_path);
    if(!metadata)
      throw std::runtime_error("Failed to open metadata output file");
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"generator\": \"multi_view_pinhole_depth_camera\",\n"
             << "  \"length_unit\": \"m\",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"object\": {\"type\": \"ellipsoid\", "
                "\"center_m\": ["
             << center.x << ", " << center.y << ", " << center.z
             << "], \"radii_m\": [" << radii.x << ", " << radii.y
             << ", " << radii.z << "]},\n"
             << "  \"camera\": {\"views\": " << options.views
             << ", \"image_width\": " << options.image_width
             << ", \"image_height\": " << options.image_height
             << ", \"vertical_fov_degrees\": "
             << options.vertical_fov_degrees << "},\n"
             << "  \"sensor_model\": {\"base_depth_noise_sigma_m\": "
             << options.base_depth_noise_sigma
             << ", \"incidence_dependent_noise\": true, "
                "\"depth_quantization_m\": "
             << options.depth_quantization
             << ", \"base_dropout_probability\": "
             << options.dropout_probability
             << ", \"grazing_angle_dropout\": true, \"outliers\": 0},\n"
             << "  \"sparsification\": {\"method\": \"voxel_reservoir\", "
                "\"voxel_size_m\": "
             << options.voxel_size << "},\n"
             << "  \"observed_noise\": {\"mean_m\": " << mean_noise
             << ", \"rms_m\": " << rms_noise << ", \"maximum_m\": "
             << maximum_noise << "},\n"
             << "  \"counts\": {\"cast_rays\": " << cast_rays
             << ", \"surface_hits\": " << surface_hits
             << ", \"accepted_measurements\": "
             << accepted_measurements << ", \"output_points\": "
             << voxels.size() << "},\n"
             << "  \"ground_truth_file\": \""
             << ground_truth_path.filename().string() << "\"\n"
             << "}\n";

    std::cout << "Simulated multi-view RGB-D point cloud\n"
              << "Length unit: m\n"
              << "Object: ellipsoid, 0.36 m x 0.26 m x 0.48 m\n"
              << "Views: " << options.views << "\n"
              << "Cast rays: " << cast_rays << "\n"
              << "Surface hits: " << surface_hits << "\n"
              << "Accepted noisy measurements: " << accepted_measurements
              << "\n"
              << "Sparse output points: " << voxels.size() << "\n"
              << "Base depth noise sigma (m): "
              << options.base_depth_noise_sigma << "\n"
              << "Observed RMS point error (m): " << rms_noise << "\n"
              << "Observed maximum point error (m): " << maximum_noise
              << "\n"
              << "Voxel size (m): " << options.voxel_size << "\n"
              << "Output: " << output.string() << "\n"
              << "Ground truth: " << ground_truth_path.string() << "\n"
              << "Metadata: " << metadata_path.string() << "\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Point-cloud simulation failed: " << error.what() << '\n';
    return 1;
  }
}
