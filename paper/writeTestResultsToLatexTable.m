function writeTestResultsToLatexTable(results, nSeg0, orders, fileName, outputDir)
% write the numbers of error and convergence rate in results to latex table

if nargin<3
    outputDir = './';
end

refRatio = 2;

sz   = size(results);
nRun = (sz(2)-1)/2+1;
N    = nSeg0*refRatio.^[0:1:nRun-1];

fileStr = [outputDir,fileName];
fid = fopen(fileStr, 'w');

format = '%6.2e';
for i=1:nRun-1
    format = [format, ' & %6.2f & %6.2e'];
end
format = [format, ' \\\\ \n'];
formSt = '%s \n';
tabFormat = 'c|c';
tabHeadLine = '$\kappa$ & ';
%Estr = '$\|\mathbf{E}_I';
Estr = '$E_{\kappa}';
for i=1:nRun-1
    tabFormat = [tabFormat,'cc'];
    tabHeadLine = [tabHeadLine, ...
        Estr,'\left(\frac{1}{', int2str(N(i)), '}\right)$ &  ${\cal O}_{\kappa}$ & '];
end
tabHeadLine = [tabHeadLine, Estr,'\left(\frac{1}{', int2str(N(nRun)), '}\right)$\\'];
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
