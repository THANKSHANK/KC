#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

// Include stb_image_write for saving PNG files.
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// -----------------------
// Basic Math Structures
// -----------------------
struct Vec2f {
    float x, y;

    // Default constructor
    Vec2f() : x(0.0f), y(0.0f) {}

    // Two-float constructor
    Vec2f(float _x, float _y) : x(_x), y(_y) {}
};

struct Vec3f {
    float x, y, z;
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float a, float b, float c) : x(a), y(b), z(c) {}
};

// -----------------------
// Utility Functions
// -----------------------
inline Vec3f normalize(const Vec3f& v) {
    float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return Vec3f(v.x / len, v.y / len, v.z / len);
}

inline float dot(const Vec3f& a, const Vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Hammersley Sequence
float RadicalInverse_VdC(unsigned int bits) {
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    return float(bits) * 2.3283064365386963e-10f; // 1/4294967296
}

Vec2f Hammersley(unsigned int i, unsigned int N) {
    return Vec2f(float(i) / float(N), RadicalInverse_VdC(i));
}

// Importance Sample GGX
Vec3f ImportanceSampleGGX(const Vec2f& Xi, const Vec3f& N, float roughness) {
    float a = roughness * roughness;
    float phi = 2.0f * 3.14159265f * Xi.x;
    float cosTheta = std::sqrt((1.0f - Xi.y) / (1.0f + (a * a - 1.0f) * Xi.y));
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);

    Vec3f H;
    H.x = sinTheta * std::cos(phi);
    H.y = sinTheta * std::sin(phi);
    H.z = cosTheta;
    // Since we assume N = (0,0,1), H is already in the same space.
    return normalize(H);
}
// Geometry Smith Function
float GeometrySmith(float roughness, float NoV, float NoL) {
    float r = roughness + 1.0f;
    float k = (r * r) / 8.0f;
    float G_V = NoV / (NoV * (1.0f - k) + k);
    float G_L = NoL / (NoL * (1.0f - k) + k);
    return G_V * G_L;
}
// BRDF Integration Function
Vec3f IntegrateBRDF(const Vec3f& V, float roughness) {
    const int sampleCount = 1024;
    float A = 0.0f;
    float B = 0.0f;
    Vec3f N(0.0f, 0.0f, 1.0f);

    for (int i = 0; i < sampleCount; i++) {
        Vec2f Xi = Hammersley(i, sampleCount);
        Vec3f H = ImportanceSampleGGX(Xi, N, roughness);
        // Reflect V about H to get L.
        float dotVH = dot(V, H);
        Vec3f L = normalize(Vec3f(
            2.0f * dotVH * H.x - V.x,
            2.0f * dotVH * H.y - V.y,
            2.0f * dotVH * H.z - V.z
        ));

        float NoL = std::max(L.z, 0.0f);
        float NoH = std::max(H.z, 0.0f);
        float VoH = std::max(dot(V, H), 0.0f);
        float NoV = std::max(dot(N, V), 0.0f);

        // Schlick Fresnel approximation.
        float Fc = std::pow(1.0f - VoH, 5.0f);
        float G = GeometrySmith(roughness, NoV, NoL);
        float G_Vis = (VoH * G) / (NoV * NoH + 1e-5f);

        // Split-sum integration.
        A += (1.0f - Fc) * G_Vis;
        B += Fc * G_Vis;
    }
    A /= sampleCount;
    B /= sampleCount;
    return Vec3f(A, B, 0.0f);
}

// -----------------------
// Main: Generate the LUT
// -----------------------
int main() {
    const int width = 256;
    const int height = 256;
    std::vector<float> lutData(width * height * 3, 0.0f);  // 3 channels: R, G, B (B remains 0)

    // Loop over a 2D grid where:
    // - X coordinate maps to NdotV in [0,1]
    // - Y coordinate maps to roughness in [0,1]
    for (int j = 0; j < height; j++) {
        float roughness = float(j) / float(height - 1);
        for (int i = 0; i < width; i++) {
            float NdotV = float(i) / float(width - 1);
            // Compute view vector V: V.z = NdotV, V.x = sqrt(1 - NdotV^2), V.y = 0
            float sinTheta = std::sqrt(std::max(1.0f - NdotV * NdotV, 0.0f));
            Vec3f V(sinTheta, 0.0f, NdotV);
            Vec3f integrated = IntegrateBRDF(V, roughness);

            int index = (j * width + i) * 3;
            lutData[index + 0] = integrated.x;  // Store term A in red
            lutData[index + 1] = integrated.y;  // Store term B in green
            lutData[index + 2] = integrated.z;  // Blue channel (unused, set to 0)
        }
    }

    // Convert the floating-point LUT data to 8-bit per channel.
    std::vector<unsigned char> outputData(width * height * 3);
    for (size_t i = 0; i < lutData.size(); i++) {
        float clamped = std::min(std::max(lutData[i], 0.0f), 1.0f);
        outputData[i] = static_cast<unsigned char>(clamped * 255.0f);
    }

    // Save the LUT as a PNG image.
    if (stbi_write_png("BRDF_LUT.png", width, height, 3, outputData.data(), width * 3)) {
        std::cout << "BRDF LUT saved as BRDF_LUT.png" << std::endl;
    }
    else {
        std::cerr << "Failed to save BRDF LUT image." << std::endl;
    }

    return 0;
}
