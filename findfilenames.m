function filenames = findfilenames(filePath,extention,specSTR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%����ĳһ���ļ��У��ҳ������׺��Ϊ".pv"���ļ����������
%filePath���ļ��е�·��
%extention���ļ���ָ����׺��   .cnt
%filecell��������������ļ���cell
%num������ͳ���ļ��ĸ�����
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