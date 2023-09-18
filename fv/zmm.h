#include "stdafx.h"

#ifndef ZMM_H
#define ZMM_H

namespace fv
{
    /// <summary>
    /// ZMM register.
    /// </summary>
    class ZMM
    {

    public:

        // Constants.

        /// <summary>
        /// Bits count.
        /// </summary>
        static const int bits = 512;

        /// <summary>
        /// Bytes count.
        /// </summary>
        static const int bytes = bits / 8;

        /// <summary>
        /// Get elements count.
        /// </summary>
        /// <typeparam name="T">Type.</typeparam>
        /// <returns>Elements count.</returns>
        template <typename T>
        static int count()
        {
            return bytes / sizeof(T);
        }

    private:

        /// <summary>
        /// Identifier.
        /// </summary>
        int id;

        /// <summary>
        /// Data.
        /// </summary>
        int8_t* data;

    public:

        // Constructors.
        
        ZMM();

        ZMM(const ZMM& z);

        ZMM& operator=(const ZMM& z);

        ZMM(ZMM&& z);

        ZMM& operator=(ZMM&& z);

        ~ZMM();

        // Access.

        /// <summary>
        /// Get id.
        /// </summary>
        /// <returns>Id.</returns>
        inline int get_id() const
        {
            return id;
        }

        /// <summary>
        /// Get element of data.
        /// </summary>
        /// <typeparam name="T">Type.</typeparam>
        /// <param name="i">Index.</param>
        /// <returns>Data element.</returns>
        template <typename T>
        T get(int i) const
        {
            const void* vp = static_cast<const void*>(&data[i * sizeof(T)]);
            const T* tp = static_cast<const T*>(vp);

            return *tp;
        }

        /// <summary>
        /// Set element of data.
        /// </summary>
        /// <typeparam name="T">Type.</typeparam>
        /// <param name="i">Index.</param>
        /// <param name="v">Value.</param>
        template <typename T>
        void set(int i,
                 T v)
        {
            void *vp = static_cast<void*>(&data[i * sizeof(T)]);
            T* tp = static_cast<T*>(vp);

            *tp = v;
        }

        /// <summary>
        /// Representation.
        /// </summary>
        /// <typeparam name="T">Type.</typeparam>
        /// <returns>Representation.</returns>
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

        /// <summary>
        /// Clear.
        /// </summary>
        void clear();
    };
}

#endif
