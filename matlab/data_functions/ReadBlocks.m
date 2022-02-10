function data = ReadBlocks(file,blocks)
temp = cell(blocks,1);
l = zeros(blocks,1);
for i=1:blocks
    temp{i} = importdata([file,'_Block',num2str(i),'.txt'],'\t');
    l(i) = length(temp{i});
end

data = cat(1,temp{:});
end

