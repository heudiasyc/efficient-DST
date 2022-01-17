#ifndef EFFICIENT_DST_MOBIUS_INVERSIONS_HPP
#define EFFICIENT_DST_MOBIUS_INVERSIONS_HPP

#include <mobius_inversion_template.hpp>


namespace efficient_DST{

	template<
		typename T,
		size_t N
	>
	struct try_linear_focal_points_computation_additive_up_inclusion : try_linear_focal_points_computation<
		T,
		N,
		up_inclusion<T, N>,
		zeta_additive_operation<T>
	>{};

	template<
		typename T,
		size_t N
	>
	struct zeta_additive_transformation_up_inclusion : zeta_transformation<
		T,
		N,
		up_inclusion<T, N>,
		zeta_additive_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_zeta_additive_up_inclusion : FMT<
		T,
		N,
		up_inclusion<T, N>,
		zeta_additive_operation<T>,
		FMT_operation_up_inclusion<T, zeta_additive_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_zeta_additive_up_inclusion : EMT<
		T,
		N,
		zeta_additive_transformation_up_inclusion<T, N>,
		up_inclusion<T, N>,
		zeta_additive_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct mobius_additive_transformation_up_inclusion : mobius_transformation<
		T,
		N,
		up_inclusion<T, N>,
		mobius_additive_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_mobius_additive_up_inclusion : FMT<
		T,
		N,
		up_inclusion<T, N>,
		mobius_additive_operation<T>,
		FMT_operation_up_inclusion<T, mobius_additive_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_mobius_additive_up_inclusion : EMT<
		T,
		N,
		mobius_additive_transformation_up_inclusion<T, N>,
		up_inclusion<T, N>,
		mobius_additive_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct try_linear_focal_points_computation_multiplicative_up_inclusion : try_linear_focal_points_computation<
		T,
		N,
		up_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>{};

	template<
		typename T,
		size_t N
	>
	struct zeta_multiplicative_transformation_up_inclusion : zeta_transformation<
		T,
		N,
		up_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_zeta_multiplicative_up_inclusion : FMT<
		T,
		N,
		up_inclusion<T, N>,
		zeta_multiplicative_operation<T>,
		FMT_operation_up_inclusion<T, zeta_multiplicative_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_zeta_multiplicative_up_inclusion : EMT<
		T,
		N,
		zeta_multiplicative_transformation_up_inclusion<T, N>,
		up_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct mobius_multiplicative_transformation_up_inclusion : mobius_transformation<
		T,
		N,
		up_inclusion<T, N>,
		mobius_multiplicative_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_mobius_multiplicative_up_inclusion : FMT<
		T,
		N,
		up_inclusion<T, N>,
		mobius_multiplicative_operation<T>,
		FMT_operation_up_inclusion<T, mobius_multiplicative_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_mobius_multiplicative_up_inclusion : EMT<
		T,
		N,
		mobius_multiplicative_transformation_up_inclusion<T, N>,
		up_inclusion<T, N>,
		mobius_multiplicative_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct try_linear_focal_points_computation_additive_down_inclusion : try_linear_focal_points_computation<
		T,
		N,
		down_inclusion<T, N>,
		zeta_additive_operation<T>
	>{};

	template<
		typename T,
		size_t N
	>
	struct zeta_additive_transformation_down_inclusion : zeta_transformation<
		T,
		N,
		down_inclusion<T, N>,
		zeta_additive_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_zeta_additive_down_inclusion : FMT<
		T,
		N,
		down_inclusion<T, N>,
		zeta_additive_operation<T>,
		FMT_operation_down_inclusion<T, zeta_additive_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_zeta_additive_down_inclusion : EMT<
		T,
		N,
		zeta_additive_transformation_down_inclusion<T, N>,
		down_inclusion<T, N>,
		zeta_additive_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct mobius_additive_transformation_down_inclusion : mobius_transformation<
		T,
		N,
		down_inclusion<T, N>,
		mobius_additive_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_mobius_additive_down_inclusion : FMT<
		T,
		N,
		down_inclusion<T, N>,
		mobius_additive_operation<T>,
		FMT_operation_down_inclusion<T, mobius_additive_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_mobius_additive_down_inclusion : EMT<
		T,
		N,
		mobius_additive_transformation_down_inclusion<T, N>,
		down_inclusion<T, N>,
		mobius_additive_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct try_linear_focal_points_computation_multiplicative_down_inclusion : try_linear_focal_points_computation<
		T,
		N,
		down_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>{};

	template<
		typename T,
		size_t N
	>
	struct zeta_multiplicative_transformation_down_inclusion : zeta_transformation<
		T,
		N,
		down_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_zeta_multiplicative_down_inclusion : FMT<
		T,
		N,
		down_inclusion<T, N>,
		zeta_multiplicative_operation<T>,
		FMT_operation_down_inclusion<T, zeta_multiplicative_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_zeta_multiplicative_down_inclusion : EMT<
		T,
		N,
		zeta_multiplicative_transformation_down_inclusion<T, N>,
		down_inclusion<T, N>,
		zeta_multiplicative_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////

	template<
		typename T,
		size_t N
	>
	struct mobius_multiplicative_transformation_down_inclusion : mobius_transformation<
		T,
		N,
		down_inclusion<T, N>,
		mobius_multiplicative_operation<T>
	>
	{};

	template<
		typename T,
		size_t N
	>
	struct FMT_mobius_multiplicative_down_inclusion : FMT<
		T,
		N,
		down_inclusion<T, N>,
		mobius_multiplicative_operation<T>,
		FMT_operation_down_inclusion<T, mobius_multiplicative_operation<T>>
	>
	{};

	template <
		typename T,
		size_t N
	>
	struct EMT_mobius_multiplicative_down_inclusion : EMT<
		T,
		N,
		mobius_multiplicative_transformation_down_inclusion<T, N>,
		down_inclusion<T, N>,
		mobius_multiplicative_operation<T>
	>
	{};

	////////////////////////////////////////////////////////////////////////////////////////////////////
}		// namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_INVERSIONS_HPP
