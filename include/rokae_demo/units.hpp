#pragma once

namespace rokae_demo
{

// Project-wide geometry convention. Input formats such as XYZ/OFF generally
// do not encode a physical unit, so callers must convert coordinates to meters
// before invoking either executable.
inline constexpr const char* kLengthUnit = "m";
inline constexpr const char* kAreaUnit = "m^2";
inline constexpr const char* kVolumeUnit = "m^3";

}  // namespace rokae_demo
