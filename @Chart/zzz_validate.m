function rPoly = validate(obj,varargin)
	% Computes validated error bounds for obj integration with a given truncation_error.
	
	% Compute F,DF, A, A_dagger, and norms
	if(nargin>1)
		trunc_error = varargin{1};
	else
		trunc_error = 0;
	end
	MN = prod(obj.modes);
	
	% finite part of F_MN(X,Y,Z)
	% add/subtract gamma has no effect on the finite part so it is omitted
	F1 = obj.var(1).diff - obj.tau*shift(obj.var(2) + obj.forcing(1));
	F2 = obj.var(2).diff - obj.tau*shift(obj.params*(obj.var(2) - obj.var(1)*obj.var(1)*obj.var(2)) - obj.var(1) + obj.forcing(2));
	F = [F1;F2];
	
	% finite parts of DF(X,Y,Z)
	% order 4 Kroneker delta. As an operator this acts as the identity.
	Kron_delta = BAoperator(eye(MN),obj.modes);
	
	% prime operator acts on vectors by Ah = h'
	Idprime = BAoperator(diag(repmat([1,1:obj.modes(1)-1],1,obj.modes(2))),obj.modes);
	
	% shift operator acts on vectors by Ah = eta(h)
	shift_subdiag = ones(1,MN);
	shift_subdiag(obj.modes(2):obj.modes(2):end) = zeros(1,obj.modes(1));
	shift_matrix = diag(shift_subdiag);
	Idshift = BAoperator([shift_matrix(end,:);shift_matrix(1:end-1,:)],obj.modes);
	
	% define finite part of DF as 4 BAoperators (components of the Jacobian)
	DF(1,1) = Idprime; %DxF1
	DF(1,2) = -obj.tau*Idshift; %DyF1
	DF(2,1) = obj.tau*Lmult(2*obj.params*shift(obj.var(1)*obj.var(2))); %DxF2
	DF(2,2) = Idprime - obj.tau*Idshift + obj.tau*Lmult(shift(obj.var(1)*obj.var(1))); %DyF2
	
	
	% compute numerical inverse of DF
	DF_matrix = [DF(1,1).matrix,DF(1,2).matrix;DF(2,1).matrix,DF(2,2).matrix];
	A_matrix = inv(DF_matrix);
	for j = 1:2
		for k = 1:2
			A(j,k) = BAoperator(A_matrix(1+(j-1)*MN:j*MN,1+(k-1)*MN:k*MN),obj.modes);
		end
	end
	
	%% Radii Polynomial Bounds
	
	%interval constants
	intval_tau = midrad(obj.tau,0);
	r = midrad(.001,0); % initial guess for Lipschitz bound
	intval_params = midrad(obj.params,0);
	a = obj.var(1).intval_norm; % interval enclosure of ||a||
	b = obj.var(2).intval_norm; % interval enclosure of ||b||
	
%             row_norms = [sum(arrayfun(@(j)A(1,j).norm,1:2));sum(arrayfun(@(j)A(2,j).norm,1:2))];
%             % non-interval computation
%             row_norms = [sum(arrayfun(@(j)A(1,j).intval_norm,1:2,'UniformOutput',false));sum(arrayfun(@(j)A(2,j).intval_norm,1:2,'UniformOutput',false))];
	row_norms(1) = A(1,1).intval_norm + A(1,2).intval_norm;
	row_norms(2) = A(2,1).intval_norm + A(2,2).intval_norm;
	normA = max(sup(row_norms));
	
	% Y_0
	y1 = A(1,1)*F(1) + A(1,2)*F(2); %[A*F]_1 %DOES THIS MULTIPLICATION HAVE TO BE IN INTERVALS OR ONLY THE NORM COMPUTATION? 
%             Y_1 = y1.norm;
	Y_1 = y1.intval_norm; 
	
	y2 = mtimes(obj.var(1),obj.var(1),Inf); %full X^2 including spillover terms
	y2 = mtimes(y2,obj.var(2),Inf); % full XXY including spillover terms
	AF_2 = A(2,1)*F(1) + A(2,2)*F(2); %convolution without spillover terms
	y2.coef(1:obj.modes(1),1:obj.modes(2)) = AF_2.coef;
%             Y_2 = y2.norm;
	Y_2 = y2.intval_norm;
	
	Y_0 = max(sup([Y_1,Y_2])) + trunc_error;
	
	% Z_0
	ADF_matrix = A_matrix*DF_matrix;
	Id_MN = eye(size(ADF_matrix));
	Id_ADF = Id_MN - ADF_matrix;
	for j = 1:2
		for k = 1:2
			I_ADF(j,k) = BAoperator(Id_ADF(1+(j-1)*MN:j*MN,1+(k-1)*MN:k*MN),obj.modes);
		end
	end
%             row_norms = [sum(arrayfun(@(j)I_ADF(1,j).norm,1:2));sum(arrayfun(@(j)I_ADF(2,j).norm,1:2))];
%             row_norms = [sum(arrayfun(@(j)I_ADF(1,j).intval_norm,1:2));sum(arrayfun(@(j)I_ADF(2,j).intval_norm,1:2))];
	row_norms(1) = I_ADF(1,1).intval_norm + I_ADF(1,2).intval_norm;
	row_norms(2) = I_ADF(2,1).intval_norm + I_ADF(2,2).intval_norm;
	Z_0 = max(sup(row_norms));
	
	% Z_1
%             Z_1 = obj.weights(1)*obj.tau/obj.modes(1)*max(1,1+obj.params*(2*obj.var(1).norm*obj.var(2).norm + obj.var(1).norm^2 + 1));
	if obj.weights == [1,1];
		Z_1 = intval_tau/obj.modes(1)*max(sup([1,1 + intval_params*(2*a*b + a^2 + 1)]));
	else
		error('not yet implemented')
	end
	
	% Z_2
%             Z_2 = 6*obj.params*r*obj.weights(1)*max(1/obj.modes(1),normA);
	Z_2 = 6*intval_params*r*max(sup([midrad(1,0)/midrad(obj.modes(1),0),normA]));
	
	% radii polynomial
	rPoly.bounds = [Z_2,Z_1,Z_0,Y_0];
	rPoly.coefs = [Z_2,Z_0 + Z_1 - 1,Y_0];
%             rPoly.radius = min(roots(rPoly.coefs));
	poly = @(x)(Z_2*x^2 + (Z_0 + Z_1 - 1)*x + Y_0);
	rPoly.radius = verifynlss(poly,.25*(1-Z_0-Z_1));
	
	if imag(rPoly.radius) ~= 0;
		sprintf('Radii polynomial has no real solutions')
	else
%                 sprintf('Maximum error: %g',rPoly.radius)
		obj.error_bound = rPoly.radius;
	end
	if nargin > 1
		precision_loss = sup(rPoly.radius)/trunc_error;
		if precision_loss > 100
			disp('Loss of precision exceeded 2 decimal places on this timestep. Consider subdividing')
		end
	end
end