#pragma once
#include <iostream>

#ifndef ENABLE_LOGGING
#define ENABLE_LOGGING 0
#endif

template <typename... Args>
constexpr void log(Args &&...args)
{
    if constexpr (ENABLE_LOGGING)
    {
        (std::cout << ... << args) << '\n';
    }
}