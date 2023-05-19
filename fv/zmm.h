#ifndef ZMM_H
#define ZMM_H

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
        int8_t data[bits];

    public:

        // Constructors.

        ZMM();

        ZMM(const ZMM& z);

        // Access.

        template <typename T>
        T get(int i) const
        {
            const void* vp = static_cast<const void*>(&data[i * sizeof(T)]);
            const T* tp = static_cast<const T*>(vp);

            return *tp;
        }

        template <typename T>
        void set(int i,
                 T v)
        {
            void *vp = static_cast<void*>(&data[i * sizeof(T)]);
            T* tp = static_cast<T*>(vp);

            *tp = v;
        }

        // Representation.

        std::string getI32Representation() const;

        std::string getF32Representation() const;

        std::string getF64Representation() const;

        // Clear.

        void clear();
    };
}

#endif
