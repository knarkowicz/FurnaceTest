#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "dds.h"

float const MATH_PI = 3.14159f;

void ExportTable(char const* name, float* lutData, unsigned lutWidth, unsigned lutHeight)
{
	FILE* f = fopen(name, "w");
	for (unsigned y = 0; y < lutHeight; ++y)
	{
		for (unsigned x = 0; x < lutWidth; ++x)
		{
			float const ndotv = (x + 0.5f) / lutWidth;
			float const roughness = (y + 0.5f) / lutHeight;
			fprintf(f, "%f, %f, %f\n", ndotv, roughness, lutData[x * 4 + y * lutWidth * 4] );
		}
	}
	fclose(f);
}

float Saturate(float x)
{
	return std::min(std::max(x, 0.0f), 1.0f);
}

float RoughnessRemap(float x)
{
	return x * x;
}

uint32_t ReverseBits(uint32_t v)
{
	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
	v = (v >> 16) | (v << 16);
	return v;
}

float SmithG1(float roughness, float ndotv)
{
	float m = RoughnessRemap( roughness );
	float m2 = m * m;
	float ndotv2 = ndotv * ndotv;
	float lambda = (-1.0f + sqrtf(1.0f + m2 * (1.0f - ndotv2) / ndotv2)) * 0.5f;
	return 1.0f / (1.0f + lambda);
}

float SmithG(float roughness, float ndotv, float ndotl)
{
	float m = RoughnessRemap( roughness );
	float m2 = m * m;
	float ndotv2 = ndotv * ndotv;
	float ndotl2 = ndotl * ndotl;
	float lambdaV = (-1.0f + sqrtf(1.0f + m2 * (1.0f - ndotv2) / ndotv2)) * 0.5f;
	float lambdaL = (-1.0f + sqrtf(1.0f + m2 * (1.0f - ndotl2) / ndotl2)) * 0.5f;
	return 1.0f / (1.0f + lambdaV + lambdaL);
}

float WeakWhiteFurnaceTest(float roughness, float ndotv)
{
	float const m = RoughnessRemap(roughness);
	float const m2 = m * m;

	float const vx = sqrtf(1.0f - ndotv * ndotv);
	float const vy = 0.0f;
	float const vz = ndotv;

	float integral = 0.0f;
	unsigned const sampleNum = 2048;
	for (unsigned i = 0; i < sampleNum; ++i)
	{
		float const e1 = (float)i / sampleNum;
		float const e2 = (float)((double)ReverseBits(i) / (double)0x100000000LL);

		float const phi = 2.0f * MATH_PI * e1;
		float const cosPhi = cosf(phi);
		float const sinPhi = sinf(phi);
		float const cosTheta = sqrtf((1.0f - e2) / (1.0f + (m2 - 1.0f) * e2)); // GGX
		float const sinTheta = sqrtf(1.0f - cosTheta * cosTheta);

		float const hx = sinTheta * cosf(phi);
		float const hy = sinTheta * sinf(phi);
		float const hz = cosTheta;

		float const vdothUnsat = vx * hx + vy * hy + vz * hz;

		float const ndoth = std::max(hz, 0.0f);
		float const vdoth = std::max(vdothUnsat, 0.0f);

		float const g1 = SmithG1(roughness, ndotv);
		float const pdf = 4.0f * vdoth / ndoth;
		integral += (g1 * pdf) / (4.0f * ndotv);
	}
	integral /= sampleNum;
	return Saturate(integral);
}

float WhiteFurnaceTest(float roughness, float ndotv)
{
	float const m = RoughnessRemap(roughness);
	float const m2 = m * m;

	float const vx = sqrtf(1.0f - ndotv * ndotv);
	float const vy = 0.0f;
	float const vz = ndotv;

	float integral = 0.0f;
	unsigned const sampleNum = 4096;
	for (unsigned i = 0; i < sampleNum; ++i)
	{
		float const e1 = (float)i / sampleNum;
		float const e2 = (float)((double)ReverseBits(i) / (double)0x100000000LL);

		float const phi = 2.0f * MATH_PI * e1;
		float const cosPhi = cosf(phi);
		float const sinPhi = sinf(phi);
		float const cosTheta = sqrtf((1.0f - e2) / (1.0f + (m2 - 1.0f) * e2)); // GGX
		float const sinTheta = sqrtf(1.0f - cosTheta * cosTheta);

		float const hx = sinTheta * cosf(phi);
		float const hy = sinTheta * sinf(phi);
		float const hz = cosTheta;

		float const vdothUnsat = vx * hx + vy * hy + vz * hz;
		float const lx = 2.0f * vdothUnsat * hx - vx;
		float const ly = 2.0f * vdothUnsat * hy - vy;
		float const lz = 2.0f * vdothUnsat * hz - vz;

		float const ndotl = std::max(lz, 0.0f);
		float const ndoth = std::max(hz, 0.0f);
		float const vdoth = std::max(vdothUnsat, 0.0f);

		float const g = SmithG(roughness, ndotv, ndotl);
		float const pdf = 4.0f * vdoth / ndoth;
		integral += (g * pdf) / (4.0f * ndotv);
	}
	integral /= sampleNum;
	return Saturate(integral);
}

int main()
{
	unsigned const LUT_WIDTH = 64;
	unsigned const LUT_HEIGHT = 64;

	float lutWeakWhiteFurnace[LUT_WIDTH * LUT_HEIGHT * 4];
	float lutWhiteFurnace[LUT_WIDTH * LUT_HEIGHT * 4];

	for (unsigned y = 0; y < LUT_HEIGHT; ++y)
	{
		float const roughness = (y + 0.5f) / LUT_HEIGHT;

		for (unsigned x = 0; x < LUT_WIDTH; ++x)
		{
			float const ndotv = (x + 0.5f) / LUT_WIDTH;

			float const whiteFurnace = WhiteFurnaceTest(roughness, ndotv);
			float const weakWhiteFurnace = WeakWhiteFurnaceTest(roughness, ndotv);

			lutWeakWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 0] = weakWhiteFurnace;
			lutWeakWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 1] = weakWhiteFurnace;
			lutWeakWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 2] = weakWhiteFurnace;
			lutWeakWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 3] = 1.0f;

			lutWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 0] = whiteFurnace;
			lutWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 1] = whiteFurnace;
			lutWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 2] = whiteFurnace;
			lutWhiteFurnace[x * 4 + y * LUT_WIDTH * 4 + 3] = 1.0f;

			printf(".");
		}
	}

	ExportTable("whiteFurnace.txt", lutWhiteFurnace, LUT_WIDTH, LUT_HEIGHT);
	ExportTable("lutWeakWhiteFurnace.txt", lutWeakWhiteFurnace, LUT_WIDTH, LUT_HEIGHT);

	SaveDDS("whiteFurnace.dds", DDS_FORMAT_R32G32B32A32_FLOAT, 16, LUT_WIDTH, LUT_HEIGHT, lutWhiteFurnace);
	SaveDDS("weakWhiteFurnace.dds", DDS_FORMAT_R32G32B32A32_FLOAT, 16, LUT_WIDTH, LUT_HEIGHT, lutWeakWhiteFurnace);
	return 0;
}