function path = getPath(input)
% path = getPath(input)
% get path of the inputs
% 
% Version: 1.0
% Date: 2016/11/24
% Author: Zhixian MA <zxma_sjtu@qq.com>

path = [];
for i = 1 : length(input)
    path = [path,'/',input{i}];
end
