#include <cstdio>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main()
{
    // Just test writing a small 2x2 image:
    unsigned char data[4 * 3] = {
        255,   0,   0,   // Red
          0, 255,   0,   // Green
          0,   0, 255,   // Blue
        255, 255,   0    // Yellow
    };

    stbi_write_png("test.png", 2, 2, 3, data, 2 * 3);
    std::printf("Wrote test.png\n");
    return 0;
}
