#ifndef FV_ZMM_H
#define FV_ZMM_H

#include "stdafx.h"

#include <string>

namespace fv
{
    // 512-bits register.
    class ZMM
    {

    public:

        // Constants.
        static const int bits = 512;
        static const int bytes = bits / 8;
        static const int i32_count = bytes / sizeof(int32_t);
        static const int f32_count = bytes / sizeof(float);
        static const int f64_count = bytes / sizeof(double);

    private:

        // Union for data access.
        union
        {
            int32_t i32[i32_count];
            float f32[f32_count];
            double f64[f64_count];
        } data;

    public:

        // Constructors.

        ZMM();

        ZMM(const ZMM& z);

        // Access.

        int32_t getI32(int i) const;

        void setI32(int i,
                    int32_t v);

        float getF32(int i) const;

        void setF32(int i,
                    float v);

        double getF64(int i) const;

        void setF64(int i,
                    double v);

        // Representation.

        std::string getI32Representation() const;

        std::string getF32Representation() const;

        std::string getF64Representation() const;

        // Clear.

        void clear();
    };
}

#endif
