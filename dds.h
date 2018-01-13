#pragma once

unsigned const DDS_FORMAT_R32G32B32A32_FLOAT = 2;

bool SaveDDS(char const* path, unsigned format, unsigned texelSizeInBytes, unsigned width, unsigned height, void const* data);