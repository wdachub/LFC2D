function writeTestResultsToLatexTableTime(results, nSeg0, orders, fileName, outputDir)
% write the numbers of error and convergence rate in results to latex table

if nargin<3
    outputDir = './';
end

refRatio = 2;

sz   = size(results);
nExact=sz(1);
nTests=sz(2);
N    = nSeg0*refRatio.^[0:1:nTests-1];

fileStr = [outputDir,fileName];
fid = fopen(fileStr, 'w');

format ='%6.2f';
for i=1:nTests-1
    format = [format, '&  %6.2f '];
end
format = [format, ' \\\\ \n'];
formSt = '%s \n';
tabFormat = 'c|c';
tabHeadLine = '$\kappa$ & ';
%Estr = '$\|\mathbf{E}_I';
Tstr = '$t_{\kappa}';
for i=1:nTests
    tabFormat = [tabFormat,'cc'];
    tabHeadLine = [tabHeadLine, ...
        Tstr,'\left(\frac{1}{', int2str(N(i)), '}\right)$ &'];
end
tabHeadLine = [tabHeadLine(1:end-1),'\\'];
fprintf(fid, formSt, '\centering');
fprintf(fid, formSt, ['\begin{tabular}{', tabFormat, '}']);
fprintf(fid, formSt, '  \hline\hline');
fprintf(fid, formSt, tabHeadLine);
fprintf(fid, formSt, '  \hline');
for i=1:sz(1)
    fprintf(fid, [int2str(orders(i)), ' & ',format], results(i,:));
end
fprintf(fid, formSt, '  \hline\hline');
fprintf(fid, formSt, '\end{tabular}');
fclose(fid);
