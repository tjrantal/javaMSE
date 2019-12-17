function indice = getIndice(array,toFind)
	indice = find(cellfun(@(x) ~isempty(strfind(lower(x),lower(toFind))),array) ==1);
	if isempty(indice)
		indice = -1;
	end
