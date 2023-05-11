#include "stdafx.h"

#include "zmm.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace fv
{
    // Constructors.

    ZMM::ZMM()
    {
        clear();
    }

    ZMM::ZMM(const ZMM& z)
    {
        for (int i = 0; i < i32_count; i++)
        {
            setI32(i, z.getI32(i));
        }
    }

    // Access.

    int32_t ZMM::getI32(int i) const
    {
        return data.i32[i];
    }

    void ZMM::setI32(int i,
                     int32_t v)
    {
        data.i32[i] = v;
    }

    float ZMM::getF32(int i) const
    {
        return data.f32[i];
    }

    void ZMM::setF32(int i,
                     float v)
    {
        data.f32[i] = v;
    }

    double ZMM::getF64(int i) const
    {
        return data.f64[i];
    }

    void ZMM::setF64(int i,
                     double v)
    {
        data.f64[i] = v;
    }

    // Representation.

    std::string ZMM::getI32Representation() const
    {
        std::stringstream ss;

        ss << "I32 [";
        ss << std::setw(10) << getI32(0);
        for (int i = 1; i < i32_count; i++)
        {
            ss << " " << std::setw(10) << getI32(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    std::string ZMM::getF32Representation() const
    {
        std::stringstream ss;

        ss << "F32 [";
        ss << std::setprecision(4) << std::setw(10) << getF32(0);
        for (int i = 1; i < f32_count; i++)
        {
            ss << " " << std::setprecision(4) << std::setw(10) << getF32(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    std::string ZMM::getF64Representation() const
    {
        std::stringstream ss;

        ss << "F64 [";
        ss << std::setprecision(4) << std::setw(10) << getF64(0);
        for (int i = 1; i < f64_count; i++)
        {
            ss << " " << std::setprecision(4) << std::setw(10) << getF64(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    // Clear.

    void ZMM::clear()
    {
        for (int i = 0; i < i32_count; i++)
        {
            setI32(i, 0);
        }
    }
}
