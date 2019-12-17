%Get matching file based on GPS file time stamp
function matchingFile = getMatchingFile(fileName, searchFolder)
	%Match is based on time stamp on the fileName 
    split = strsplit(fileName,'.txt');
    usplit = strsplit(split{1},'_');
    
	dateTimeToMatch = [usplit{end-1} '_' usplit{end}(1:end-2)]; %Date string including minutes
	dNum = datenum(dateTimeToMatch,'yyyy-mm-dd_HHMM');
	minInc = 1/(24*60);
	%Create time stamp strings that are off by -2 to 2 min
	matchDates = cellfun(@(x) datestr(x,'yyyy-mm-dd_HHMM'),num2cell(dNum+[-2*minInc:minInc:2*minInc]),'uniformoutput',false);
	%Look through match time stamps in 0 +1 -1 ... order
	matchDates = matchDates([3,4,2,5,1]);
	%Get list of files to match with respect to the time stamp
	fList = dir([searchFolder '/*.txt']);
	matchingFileTest = cellfun(@(x) getIndice({fList(:).name},x),matchDates);
	matchingFileIndex = matchingFileTest(find(matchingFileTest > 0,1,'first'));
	%Return a single matching file name
	matchingFile = fList(matchingFileIndex).name;
	
