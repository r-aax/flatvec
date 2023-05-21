#include "stdafx.h"

#ifndef ZMM_H
#define ZMM_H

namespace fv
{
    // 512-bits register.
    class ZMM
    {

    public:

        // Constants.

        static const int bits = 512;
        static const int bytes = bits / 8;

        // Count of elements.
        template <typename T>
        static int count()
        {
            return bytes / sizeof(T);
        }

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
        template <typename T>
        std::string get_representation() const
        {
            std::stringstream ss;

            ss << "[";
            ss << std::setw(10) << get<T>(0);
            for (int i = 1; i < count<T>(); i++)
            {
                ss << " " << std::setw(10) << get<T>(i);
            }
            ss << "]";

            return std::string {ss.str()};
        }

        // Clear.
        void clear();
    };
}

#endif
