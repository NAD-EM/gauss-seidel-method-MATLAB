%Calculates solution for a linear system of equations using the Gauss–Seidel method:
%matrix A must be diagonally dominant.
%Usage example: "s = GaussSeidel([8 1 6 ; 4 9 2 ; 3 5 7], [0.5 0.2 0.3], 0.00001)"

%{
	GaussSeidel: script to calculate a solution for a linear system of equations using the Gauss–Seidel method.
    Copyright (C) 2017  NAD-EM

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

    You can contact the author through email:
	<git.nad.em.00@gmail.com>.
	Or through their GitHub profile:
	<https://github.com/NAD-EM>.
%}

function s = GaussSeidel(A, res, TolerableErrorRange)
format long;

[m,n] = size(A); % rows / columns
s = zeros(1, n); % solution to the system

temp = 0;
CurrentError = 100;

while CurrentError > TolerableErrorRange
    prev = s; % save previous aproximation
    columns = 0;
    for rows = 1 : m % iterate thru the rows
        columns = columns + 1;
        for col = 1 : n % iterate coefficients (columns) of current row
            if col ~= columns
                temp = temp - ( A(rows, col) * s(col) );
            else
                temp = temp + res(col);
            end
        end
        s(columns) = temp / A(columns, columns);
        temp = 0;
    end
    CurrentError = abs( (s(n) - prev(n)) / s(n) ) * 100;
end

end
