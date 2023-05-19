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
            set<int32_t>(i, z.get<int32_t>(i));
        }
    }

    // Representation.

    std::string ZMM::getI32Representation() const
    {
        std::stringstream ss;

        ss << "I32 [";
        ss << std::setw(10) << get<int32_t>(0);
        for (int i = 1; i < i32_count; i++)
        {
            ss << " " << std::setw(10) << get<int32_t>(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    std::string ZMM::getF32Representation() const
    {
        std::stringstream ss;

        ss << "F32 [";
        ss << std::setprecision(4) << std::setw(10) << get<float>(0);
        for (int i = 1; i < f32_count; i++)
        {
            ss << " " << std::setprecision(4) << std::setw(10) << get<float>(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    std::string ZMM::getF64Representation() const
    {
        std::stringstream ss;

        ss << "F64 [";
        ss << std::setprecision(4) << std::setw(10) << get<double>(0);
        for (int i = 1; i < f64_count; i++)
        {
            ss << " " << std::setprecision(4) << std::setw(10) << get<double>(i);
        }
        ss << "]";

        return std::string {ss.str()};
    }

    // Clear.

    void ZMM::clear()
    {
        for (int i = 0; i < i32_count; i++)
        {
            set<int32_t>(i, 0);
        }
    }
}
