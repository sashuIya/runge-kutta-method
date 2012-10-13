Overview
--------

This program solves the differential equation system 

	x' = x + y
	y' = (alpha - 1) * y + (alpha - 2) * x - x^3 - x^2 * y

with Runge-Kutta method.

Usage:
------

1. Compile with g++
2. Input alpha, x0, y0 and tau (the size of Runge-Kutta step)

Notice, that it makes steps only in one direction (depends on sign(tau)).

License
-------
    Copyright (C) 2012  Alexander Lapin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
Contacts:
---------
Alexander Lapin, lapinra@gmail.com
