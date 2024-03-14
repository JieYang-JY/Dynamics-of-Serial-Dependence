function [p_FM, chisquare] = getFishersP(plist)

adj_plist = plist .* 2;
temp_p_list = []; 
for i = 1:length(adj_plist); temp_p_list = [temp_p_list -2*log(adj_plist(i))]; end
df = 2 * length(plist); chisquare = sum(temp_p_list);
format shortE; p_FM = 1 - chi2cdf(chisquare, df);

