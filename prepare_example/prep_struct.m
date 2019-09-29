pt = rmfield(pt,'ieeg_name');
pt = rmfield(pt,'dmin');
pt = rmfield(pt,'elecNotes');
pt = rmfield(pt,'chLocationFile');
pt = rmfield(pt,'dmin');
pt = rmfield(pt,'electrode_labels');
pt = rmfield(pt,'stats');
pt = rmfield(pt,'data');
pt = rmfield(pt,'seq_matrix');
pt = rmfield(pt,'removed');
pt = rmfield(pt,'newSOZNames');
pt = rmfield(pt,'resecElecs');
pt = rmfield(pt,'resecLabels');

for i = 1:length(pt)
    pt(i).allTimes = [];
    pt(i).runTimes = [];
    pt(i).chunkFiles = {};
    
end