// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

#include "precomp.h"
#include "inc/ColorFix.hpp"

#include <DirectXPackedVector.h>

static DirectX::XMVECTOR XM_CALLCONV cbrtf_est_simd(DirectX::FXMVECTOR a)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
        cbrtf(a.vector4_f32[0]),
        cbrtf(a.vector4_f32[1]),
        cbrtf(a.vector4_f32[2]),
        cbrtf(a.vector4_f32[3]),
    } } };
    return Result.v;
#else

    // http://metamerist.com/cbrt/cbrt.htm showed a great estimator for the cube root:
    //   float_as_uint32_t / 3 + 709921077
    // It's similar to the well known "fast inverse square root" trick. Lots of numbers around 709921077 perform
    // at least equally well to 709921077, and it is unknown how and why 709921077 was chosen specifically.
    //
    // Integer division is expensive, so we do it the compiler way: Multiplicate and shift.
    //   x / 3 == (x * 0xaaaaaaab) >> 1
    //
    // Unfortunately DirectXMath lacks integer routines so we have to do that ourselves.
#if defined(_XM_ARM_NEON_INTRINSICS_)
    const auto div3 = vshrq_n_u32(vmulq_n_u32(vreinterpretq_u32_f32(a), 0xaaaaaaab), 1);
    // NEON is so confusing... Its fairly sparse documentation lists a vaddq_n_u32, but it doesn't seem to exist anywhere.
    const auto x = vreinterpretq_f32_u32(vaddq_u32(div3, vdupq_n_u32(709921077)));
#else
    const auto div3 = _mm_srli_epi32(_mm_mul_epu32(_mm_castps_si128(a), _mm_set1_epi32(0xaaaaaaab)), 1);
    const auto x = _mm_castsi128_ps(_mm_add_epi32(div3, _mm_set1_epi32(709921077)));
#endif

    // One round of Halley's method. It follows the Wikipedia article at
    //   https://en.wikipedia.org/wiki/Cube_root#Numerical_methods
    // For `a`s in the range between 0 and 1, this has a maximum error of 1.6e-5.
    // Newton's method would've been sufficient for our needs as well, but strangely enough it
    // ran a tiny bit slower in practice when compiled with MSVC. Given that the difference
    // was <1ns I didn't try to figure out why and just used the more precise method.
    const auto x3 = DirectX::XMVectorMultiply(DirectX::XMVectorMultiply(x, x), x);
    const auto dividend = DirectX::XMVectorMultiply(x, DirectX::XMVectorAdd(x3, DirectX::XMVectorAdd(a, a)));
    const auto divisor = DirectX::XMVectorAdd(DirectX::XMVectorAdd(x3, x3), a);
    return DirectX::XMVectorDivide(dividend, divisor);
#endif
}

static DirectX::XMVECTOR XM_CALLCONV linearSrgbToOklab(DirectX::PackedVector::XMCOLOR color)
{
    static const DirectX::XMMATRIX lmsM{
        { 0.4122214708f, 0.2119034982f, 0.0883024619f, 0.0f },
        { 0.5363325363f, 0.6806995451f, 0.2817188376f, 0.0f },
        { 0.0514459929f, 0.1073969566f, 0.6299787005f, 0.0f },
        { 0.0f, 0.0f, 0.0f, 0.0f },
    };
    static const DirectX::XMMATRIX labM{
        { 0.2104542553f, 1.9779984951f, 0.0259040371f, 0.0f },
        { 0.7936177850f, -2.4285922050f, 0.7827717662f, 0.0f },
        { -0.0040720468f, 0.4505937099f, -0.8086757660f, 0.0f },
        { 0.0f, 0.0f, 0.0f, 0.0f },
    };

    const auto c = DirectX::PackedVector::XMLoadColor(&color);
    const auto lms = DirectX::XMVector3Transform(c, lmsM);
    const auto lms_ = cbrtf_est_simd(lms);
    return DirectX::XMVector3Transform(lms_, labM);
}

static DirectX::PackedVector::XMCOLOR XM_CALLCONV oklabToLinearSrgb(DirectX::XMVECTOR lab)
{
    static const DirectX::XMMATRIX lmsM{
        { 4.0767416621f, -1.2684380046f, -0.0041960863f, 0.0f },
        { -3.3077115913f, 2.6097574011f, -0.7034186147f, 0.0f },
        { 0.2309699292f, -0.3413193965f, 1.7076147010f, 0.0f },
        { 0.0f, 0.0f, 0.0f, 0.0f },
    };
    static const DirectX::XMMATRIX labM{
        { 1.0f, 1.0f, 1.0f, 0.0f },
        { 0.3963377774f, -0.1055613458f, -0.0894841775f, 0.0f },
        { 0.2158037573f, -0.0638541728f, -1.2914855480f, 0.0f },
        { 0.0f, 0.0f, 0.0f, 0.0f },
    };

    const auto lms_ = DirectX::XMVector3Transform(lab, labM);
    const auto lms = DirectX::XMVectorMultiply(DirectX::XMVectorMultiply(lms_, lms_), lms_);
    const auto c = DirectX::XMVector3Transform(lms, lmsM);
    DirectX::PackedVector::XMCOLOR color;
    DirectX::PackedVector::XMStoreColor(&color, c);
    return color;
}

COLORREF GetPerceivableColor(COLORREF reference, COLORREF color) noexcept
{
    const auto oklabReference = linearSrgbToOklab(reference);
    const auto oklabColor = linearSrgbToOklab(color);
    const auto deltas = DirectX::XMVectorSubtract(oklabReference, oklabColor);
    const auto delta = DirectX::XMVector3LengthSq(deltas);

    if (DirectX::XMVectorGetX(delta) >= 0.04f)
    {
        return color;
    }

    // TODO: broken
    auto newColor = oklabColor;
    if (oklabReference.m128_f32[0] >= 0.5f)
    {
        newColor.m128_f32[0] = oklabReference.m128_f32[0] - sqrtf(-25.0f * deltas.m128_f32[1] * deltas.m128_f32[1] - 25.0f * deltas.m128_f32[2] * deltas.m128_f32[2] + 1) * 0.2f;
    }
    else
    {
        newColor.m128_f32[0] = oklabReference.m128_f32[0] + sqrtf(-25.0f * deltas.m128_f32[1] * deltas.m128_f32[1] - 25.0f * deltas.m128_f32[2] * deltas.m128_f32[2] + 1) * 0.2f;
    }

    return oklabToLinearSrgb(newColor) | (color & 0xff000000);
}
