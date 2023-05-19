#ifndef ARRAY_MANAGER_H
#define ARRAY_MANAGER_H

#include "stdafx.h"

namespace fv
{
    template<typename T>
    class ArrayManager
    {

    private:

        T* data;
        int size;

    public:

        explicit ArrayManager<T>(int s)
            : size {s}
        {
            data = new T[size];
        }

        ~ArrayManager<T>()
        {
            delete[] data;
        }

        T* get_data()
        {
            return data;
        }

        void print() const
        {
            std::cout << "Array (" << size << "):";

            for (int i = 0; i < size; i++)
            {
                std::cout << " " << data[i];
            }

            std::cout << std::endl;
        }

        void generate_random(T lo,
                             T hi);

        T maxDiff(ArrayManager<T>& am) const
        {
            T md {};

            for (int i = 0; i < size; i++)
            {
                T d = abs(data[i] - am.get_data()[i]);

                md = std::max(md, d);
            }

            return md;
        }
    };
}

#endif
