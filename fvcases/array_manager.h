#include "stdafx.h"

#ifndef ARRAY_MANAGER_H
#define ARRAY_MANAGER_H

namespace fv
{
    /// <summary>
    /// Array manager.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    template<typename T>
    class ArrayManager
    {

    private:

        /// <summary>
        /// Data.
        /// </summary>
        T* data;

        /// <summary>
        /// Size.
        /// </summary>
        int size;

    public:

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="s">Size.</param>
        explicit ArrayManager<T>(int s)
            : size {s}
        {
            data = new T[size];

            for (int i = 0; i < s; i++)
            {
                data[i] = static_cast<T>(0);
            }
        }

        /// <summary>
        /// Constructor from file name.
        /// </summary>
        /// <param name="s">Size.</param>
        /// <param name="fn">File name.</param>
        /// <param name="skip_one_value">Skip one value between reads.</param>
        explicit ArrayManager<T>(int s,
                                 std::string fn,
                                 bool skip_one_value=false)
            : ArrayManager<T>(s)
        {
            std::ifstream f(fn);

            if (f)
            {
                char c;

                for (int i = 0; i < s; i++)
                {
                    if (!f.eof())
                    {
                        f >> data[i] >> c;

                        if (skip_one_value)
                        {
                            if (!f.eof())
                            {
                                T v;

                                f >> v >> c;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        break;
                    }
                }

                f.close();
            }
        }

        /// <summary>
        /// Destructor.
        /// </summary>
        ~ArrayManager<T>()
        {
            delete[] data;
        }

        /// <summary>
        /// Get array size.
        /// </summary>
        /// <returns>Size.</returns>
        int get_size() const
        {
            return size;
        }

        /// <summary>
        /// Get data.
        /// </summary>
        /// <returns>Data.</returns>
        T* get_data()
        {
            return data;
        }

        /// <summary>
        /// Print.
        /// </summary>
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

        /// <summary>
        /// Maximum diff with another array.
        /// </summary>
        /// <param name="am">Another array.</param>
        /// <returns>Max difference.</returns>
        T max_diff(ArrayManager<T>& am) const
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
