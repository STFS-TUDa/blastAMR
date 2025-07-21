/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef DDD_H
#define DDD_H
#include "demandDrivenData.H"
#include <memory>
#include <type_traits>

// Pointer resource management in a backwards compatible manner
template <typename T>
struct HasResetWithNullPtr {
private:
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().reset(nullptr), std::true_type{});

    template <typename>
    static std::false_type test(...);
public:
    static constexpr bool value = decltype(test<T>(0))::value;
};
template<typename T>
void deleteDemandDrivenDataOrResetImpl(T& t, std::true_type) {
    t.reset(nullptr);
}
template<typename T>
void deleteDemandDrivenDataOrResetImpl(T& t, std::false_type) {
    deleteDemandDrivenData(t);
}
template<typename T>
void deleteDemandDrivenDataOrReset(T& t) {
    deleteDemandDrivenDataOrResetImpl(t, std::integral_constant<bool, HasResetWithNullPtr<T>::value>{});
}

// For future reference, this can be simplified
//template<class T>
//void deleteDemandDrivenDataOrReset(T& t) {
//    if constexpr (requires {t.reset(nullptr);})
//    {
//        t.reset(nullptr);
//    } else {
//        deleteDemandDrivenData(t);
//    }
//}

// Assigning raw and smart pointers; managing ownership in a backwards compatible manner
template<typename Target, typename Source>
struct isPtrAssignableHelper {
private:
    template<typename T, typename S>
    static auto test(int) -> decltype(std::declval<T&>() = std::declval<S>(), std::true_type{});

    template<typename, typename>
    static std::false_type test(...);
public:
    static constexpr bool value = decltype(test<Target, Source>(0))::value;
};
template<typename TargetPtr, typename SourcePtr>
void assignPtrCompatImpl(TargetPtr& target, SourcePtr& source, std::true_type) {
    target = source;
    source = nullptr;
}
template<typename TargetPtr, typename SourcePtr>
void assignPtrCompatImpl(TargetPtr& target, SourcePtr& source, std::false_type) {
    target = source.release();
}
template<typename TargetPtr, typename SourcePtr>
void assignPtrCompat(TargetPtr& target, SourcePtr& source) {
    assignPtrCompatImpl(target, source, std::integral_constant<bool, isPtrAssignableHelper<TargetPtr, SourcePtr>::value>{});
}

template<typename TargetPtr, typename SourcePtr>
void assignUniquePtrCompatImpl(TargetPtr& target, SourcePtr& source, std::true_type) {
    target = source;
    source = nullptr;
}

// When not assignable
template<typename TargetPtr, typename SourcePtr>
void assignUniquePtrCompatImpl(TargetPtr& target, SourcePtr& source, std::false_type) {
    target.reset(source);
    source = nullptr;
}

template<typename TargetPtr, typename SourcePtr>
void assignUniquePtrCompat(TargetPtr& target, SourcePtr& source) {
    assignUniquePtrCompatImpl(target, source, std::integral_constant<bool, isPtrAssignableHelper<TargetPtr, SourcePtr>::value>{});
}

// For future reference, this can also be simplified
//template<typename TargetPtr, typename SourcePtr>
//void assignPtrCompat(TargetPtr& target, SourcePtr& source)
//{
//    if constexpr (
//        std::is_assignable_v<
//                std::remove_cvref_t<decltype(target)>,
//                std::remove_cvref_t<decltype(source)>>
//    ) {
//        target = source;
//        source = nullptr;
//    }
//    else
//    {
//        target = source.release();
//    }
//}
//template<typename TargetPtr, typename SourcePtr>
//void assignUniquePtrCompat(TargetPtr& target, SourcePtr& source)
//{
//    if constexpr (
//        std::is_assignable_v<
//                std::remove_cvref_t<decltype(target)>,
//                std::remove_cvref_t<decltype(source)>>
//    ) {
//        target = source;
//    }
//    else
//    {
//        target.reset(source);
//    }
//    source = nullptr;
//}

#endif
