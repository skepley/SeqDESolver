classdef Scalar
%SCALAR - SeqDE class for representing analytic scalar functions.
%
% The Scalar class is a finite approximation for an analytic scalar of the form, f: D --> R^n where D is
% a d-dimensional polydisc in C^d. Scalar representations are specified as double or intval coefficients
% with respect to a Taylor, Fourier, or Chebyshev series expansion.
%
%   SCALAR constructor syntax:
%       ScalarObj = Scalar()
%
%   SCALAR properties:
%       Property 1 - description
%       Property 2 - description
%
%   SCALAR methods:
%       Method 1 - description
%       Method 2 - description
%
%   Examples:
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Subclasses: none
%   Superclasses: none
%   Other classes required: none
%   Other m-files required: INTLAB toolbox
%   MAT-files required: none
%
%   See also: @BAscalar (previous version of this class), @Chart

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Jun-2016; Last revision: 08-Aug-2018
%
%   Revision History:
%   11-Jul-2017 - Support for interval and Scalar coefficients added
%   15-Aug-2017 - Reverted FFT based convolution to a classical algorithm. FFT is faster but numerically unstable, especially for intval coefficients
%   08-Aug-2018 - Class completely overhauled and renamed from BAscalar to SCalar. Numerous improvements including:
%       full class code refactorization and organization of class folder
%       class inheritance changed from handle to value
%       support for additional bases (Fourier or Chebyshev)
%       Intlab polynom dependency removed
%       support for higher dimensions
%       streamlined truncation methods
%       subscript reference methods
%       numerous efficiency gains, bug fixes, and improvements in data type consistency
%
%   ToDo:
%   Merge and simplify evaluation, fixtime, derivative, and truncation methods/function calls
%   Rigorous evaluation of sin/cos/exp, inverses, and fractional powers
%   Finish vectorization of all methods
%   Finish high precision inteval FFT and reimplement FFT-based fast convolution
%   Add polynomial evaluation over the Banach algebra
    
    %% Properties
    properties
    end % end properties
    
    %% Methods
    methods
        function obj = Scalar(input1,input2,input3,varargin)
            %SCALAR - class constructor
            
        end % end class constructor
    end % end methods
    
end % end class
