function [cartesianComponent,varargout] = complx2cart(input)
    % complex2cart - decompose the complex input signal into real and imag part
    %
    % Syntax:  [cartesianComponent] = complex2cart(input)
    %
    % Inputs:
    %    input : a column vector [Nx1]
    %
    % Outputs:
    %    cartesianComponent : the real and imag decompostion concatenated into a single 2-column matrix
    %                         if the function is called with a second output argument, columns are separated
    %                         into two different column vectors
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: complex()
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020
    %------------- BEGIN CODE --------------
    
    sizeIn = size(input);
    if prod(sizeIn(1:2))~=sizeIn(1)
        error('TELECOM205:complex2cart:inputType','Argument "input" should be a (single) column vector')
    end
    
    if nargout == 1 
        cartesianComponent = [real(input) imag(input)];
    elseif nargout == 2
        cartesianComponent  = real(input);
        varargout{1}        = imag(input);
    else
        error('TELECOM205:complex2cart:outputNargout','Number of output arguments must be either one or two')
    end





    %------------- END OF CODE --------------