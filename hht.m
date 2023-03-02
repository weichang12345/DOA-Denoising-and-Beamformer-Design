function y = hht(x, t, thr)
% hht is the result of Hilbert-Huang Transform. Where x is the input signal, t is timeline, thr is the threshold when
% checking the function is an IMF or not at step 7, y is the output, where
% each row of y is one of the IMFs of x.

y = [];
%% Step 1 : Initial condition
input = x;
x0 = x;  %% Because c is 0 in the beginning
k = 0;   %% Limit counter, the value can be set arbitrary in order to prevent endless loop

while 1
    while 1
        %% Step 2 : Find the local peaks
        peakx_input = [t(1), t(islocalmax(input)), t(end)];
        peaky_input = [input(1), input(islocalmax(input)), input(end)];
        %% Step 3 : Connect local peaks
        emax = spline(peakx_input, peaky_input, t);
        %% Step 4 : Find the local dips
        dipx_input = [t(1), t(islocalmin(input)), t(end)];
        dipy_input = [input(1), input(islocalmin(input)), input(end)];
        %% Step 5 : Connect local dips
        emin = spline(dipx_input, dipy_input, t);
        %% Step 6-1 : Compute the mean
        z = (emax + emin)/2;
        %% Step 6-2 : Compute the residue
        h = input - z;
        %% Step 7 : Check whether h is an IMF
        peakx_h = t(islocalmax(h));
        peaky_h = h(islocalmax(h));
        dipx_h = t(islocalmin(h));
        dipy_h = h(islocalmin(h));
        if(numel(dipx_h)<=1 || numel(dipy_h)<=1)  % If h has less than 2 extremas, skip the interpolation and redo the loop
            input = h;
            k = k+1;
            break
        end
        u1 = spline(peakx_h, peaky_h, t);
        u0 = spline(dipx_h, dipy_h, t);
        if(min(peaky_h) > 0 && max(dipy_h) < 0 && all(abs((u1+u0)/2) < thr) || k == 1000)  % k can be set arbitrary in order to prevent endless loop
            y = [y; h];
            %% Step 8-1 : Calculate x0
            x0 = x0 - h;
            break
        else
            input = h;
            k = k+1;
        end
    end
    %% Step 8-2 : Check whether x0 is a function with no more than 3 extreme points
    if(numel(x0(islocalmax(x0))) + numel(x0(islocalmin(x0))) < 4 || k == 1000)  % k can be set arbitrary in order to prevent endless loop
        break
    else
        input = x0;
    end
end

y = [y; x0];

end