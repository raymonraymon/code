// ***********************************************************************
// Assembly         : SpaceFrame
// Author           : chenruishan
// Created          : 09-28-2014
//
// Last Modified By : chenruishan
// Last Modified On : 09-19-2014
// ***********************************************************************
// <copyright file="Displace.h" company="">
//     Copyright (c) . All rights reserved.
// </copyright>
// <summary></summary>
// ***********************************************************************
#pragma once

/// <summary>
/// Class CDisplace.
/// </summary>
class CDisplace
{
public:
	CDisplace(void);
	CDisplace(const int &No, double* array, bool* arraybool);
	~CDisplace(void);
	void show(void);
public:
	/// <summary>
	/// The m_ no
	/// </summary>
	int	   m_No;
	/// <summary>
	/// The m_array
	/// </summary>
	double m_array[6];
	/// <summary>
	/// The m_arraybool
	/// </summary>
	bool   m_arraybool[6];
};
