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

        T* getData()
        {
            return data;
        }

        T maxDiff(ArrayManager<T>& am) const
        {
            T md {};

            for (int i = 0; i < size; i++)
            {
                T d = abs(data[i] - am.getData()[i]);

                md = std::max(md, d);
            }

            return md;
        }
    };
}

#endif
