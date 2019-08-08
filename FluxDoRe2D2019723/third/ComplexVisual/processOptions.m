function [varargout]=processOptions (opt,names,stdvalues)
% Optionen bzw. Default-Value aus struct lesen
% Durchläuft die vorgegebenen Feldnamen in "names" und
% gibt (in dieser Reihenfolge) den Standardwert aus "stdvalues" 
% zurück, sofern in der Structure "opt" nicht ein 
% anderer Wert steht.
% Input:
%     opt         Structure mit Optionen
%     names       Cell-Array mit zu suchenden Feldern 
%     stdvalues   Cell-Array mit Standardbelegung für "names"
% Output:
%     variabel    liefert genausoviele Ausgabewerte, wie Felder in
%                 Cell-Array "names", wobei der Wert vom jeweiligen
%                 "stdvalue" verwendet wird, sofern in "opt" nicht
%                 ein anderer Wert spezifiziert worden ist.
%
% See also 
%   COMPLEXVISUALLICENSE
%

varargout=cell(1,length(names));

for i=1:length(names)
   if (isfield(opt,names{i}))
      varargout{i}=opt.(names{i}); 
   else
      varargout{i}=stdvalues{i};
   end
end
end