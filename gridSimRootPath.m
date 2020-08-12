function rootPath = gridSimRootPath
	[pathDir] = fileparts(mfilename('fullpath'));
	if contains(mfilename('fullpath'), 'LiveEditorEvaluation', 'IgnoreCase', true)
		[pathDir] = fileparts(matlab.desktop.editor.getActiveFilename);
	end
	if isempty(pathDir)
		[pathDir] = fileparts(which(taskName));
		if isempty(pathDir),    pathDir     = pwd;  end
	end
	rootPath = pathDir;
end