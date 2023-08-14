function filenames = findfilenames(filePath,extention,specSTR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%遍历某一个文件夹，找出里面后缀名为".pv"的文件，并输出。
%filePath是文件夹的路径
%extention是文件的指定后缀名   .cnt
%filecell是用来保存输出文件的cell
%num是用来统计文件的个数。
%author: 476443785@qq.com

% if nargin<3;
%     specSTR=[]; 
% end

% filePath='E:\Users\Jack\Desktop\my_study\anjiabao_20121207'
% extention='.cnt'

files=dir(filePath);
fileLen=length(files);
filenames{1}=' ';
num=1; 

if nargin==3;
for i=1:fileLen    
    if files(i).isdir~=1 && ~isempty(strfind(files(i).name,extention)) && ~isempty(strfind(files(i).name,specSTR))
        %files(i).name
        filenames{num}=files(i).name;   
        num=num+1;
    end
end

elseif nargin==2
for i=1:fileLen    
    if files(i).isdir~=1 && ~isempty(strfind(files(i).name,extention))
        %files(i).name
        filenames{num}=files(i).name;
        num=num+1;
    end
end

else
    error('input arguements error')
end

end