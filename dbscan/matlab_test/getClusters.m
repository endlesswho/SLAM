function [cluster, centroids] = getClusters( points, IDX )
%GETCLUSTER ���� dbscan �ľ��������������ݺ�ÿ����ص�����
%   {return} cluster: 1*n��cell���͵ĵ������,ÿ��Ԫ��Ϊһ����ص����е㼯
%            centroids: ��Ӧÿ����ص���������
pointsNum = size(IDX, 1);
clusterNum = size(unique(IDX), 1) - 1;
cluster_ = cell(1, clusterNum);

for i = 1:pointsNum
    idx = IDX(i);
    if idx ~= 0
    cluster_{idx} = [cluster_{idx} points(i,:)'];
    end
end
%% filter out small cluster
% cluster = cluster_;
cluster = cell(1); % expand automatically
ii = 1;
for j = 1:clusterNum
    if length(cluster_{j}) > 100
        cluster{ii} = cluster_{j};
        ii = ii + 1;
    end
end
%% calculate centroids
validClusterNum = length(cluster);
centroids = zeros(3, validClusterNum);
for k = 1:validClusterNum
    centroids(:,k) = mean(cluster{k}, 2);
end
end

