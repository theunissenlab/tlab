function [p_forward, dwt, dwf, best_dwt, best_dwf, cm_dwt, cm_dwf, quad_assym, sep_q1, sep_q2]= filter_transfer(forward, samprate, fstep);

% Get dimensions
[nb nt] = size(forward);

% Smoothing filter
s_forward = fliplr(forward);  % The filter must be flipped to get the transfer function

% Window in time and frequency.
%wht = hanning(nt);
%for i=1:nb
    %s_forward(i,:) = s_forward(i,:).*wht';
    %end
%whf = hanning(nb);
%for i=1:nt
%    s_forward(:,i) = s_forward(:,i).*whf;
%end

sf_forward = fftshift(fft2(s_forward));
p_forward = sf_forward.*conj(sf_forward);

% Obtain labels for time and frequency axis
fcircle = 1/fstep;
for i=1:nb
    dwf(i) = ((i-1)/nb)*fcircle;
    if ( dwf(i) > fcircle/2 )
        dwf(i) = dwf(i) - fcircle;
    end
end
dwf = fftshift(dwf);
if ( dwf(1) > 0 ) 
    dwf(1) = -dwf(1);
end


fcircle = samprate;
for i=1:nt
    dwt(i) = ((i-1)/nt)*fcircle;
    if ( dwt(i) > fcircle/2 )
        dwt(i) = dwt(i) - fcircle;
    end
end
dwt = fftshift(dwt);
if ( dwt(1) > 0 )
    dwt(1) = -dwt(1);
end




% Find best modulation frequency from filter
max_forward = max(max(p_forward));
[ind_dwf ind_dwt] = find(p_forward == max_forward);

best_dwt = dwt(ind_dwt(1));
best_dwf = dwf(ind_dwf(1));
if (best_dwt*best_dwf > 0 )
    best_dwt = abs(best_dwt);
    best_dwf = abs(best_dwf);
else
    best_dwf = abs(best_dwf);
    best_dwt = -abs(best_dwt);
end

% Caculate relative power in first and second quadrant and center of mass
% for dwt and dwf
dwt0 = find(dwt == 0 );
dwf0 = find(dwf == 0 );
q1 = p_forward(dwf0:end, dwt0:end);
q2 = p_forward(dwf0:end, 1:dwt0 );
cm_dwf = dwf(dwf0:end)*sum(q1,2) + dwf(dwf0:end)*sum(q2,2);
cm_dwt = dwt(dwt0:end)*sum(q1,1)' + abs(dwt(1:dwt0))*sum(q2,1)';

p1 = sum(sum(q1));
p2 = sum(sum(q2));
p0 = sum(q2(:,1));

cm_dwf = cm_dwf./(p1+p2-p0);
cm_dwt = cm_dwt./(p1+p2-p0);

quad_assym = (p2-p1)/(p2+p1);

[u, s, v] = svd(q1);
s_diag = diag(s);
sep_q1 = s_diag(1).^2/sum(s_diag.^2);

[u, s, v] = svd(q2);
s_diag = diag(s);
sep_q2 = s_diag(1).^2/sum(s_diag.^2);
