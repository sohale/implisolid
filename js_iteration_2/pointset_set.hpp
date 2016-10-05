/*
    Used for debugging purposes and visualisation for sake of debugging.
*/
#pragma once

constexpr bool STORE_POINTSETS = true;

std::map< std::string, vectorized_vect > point_set_set;

// LOG_POINTSET
inline void STORE_POINTSET(std::string key_name, const vectorized_vect & centroids)
{
    if (STORE_POINTSETS) {
        vectorized_vect ps2 = centroids;  // make a copy
        if (VERBOSE_QEM) {
            clog << "STORE_POINTSET" << key_name << ": "  // "3.centroids: "
                << ps2.shape()[0] << "x" << ps2.shape()[1] << " : " << ps2[0][0] << std::endl;
        }
        auto number_erased = point_set_set.erase(std::string(key_name));
        if (VERBOSE_QEM) {
            clog << "number_erased : " << number_erased << std::endl;
        }
        point_set_set.emplace(std::make_pair(std::string(key_name), ps2));  // uses move constructor
    }
}
