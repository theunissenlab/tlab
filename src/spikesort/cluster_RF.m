function labels = cluster_RF(data, k, handles)

if get(handles.spectralcluster, 'Value') == 1

    W = SimGraph(data, 1);
    A = SpectralClustering(W, k, 3);
    labels = full(sum(bsxfun(@times, A, 1:6),2));
elseif get(handles.kmeanscluster, 'Value') == 1
    
    if size(data,1)<size(data,2)
        k=min(k, size(data,1));
        labels = fkmeans(data,k);
    else
        labels = kmeans(data,k, 'emptyaction', 'singleton');
    end
elseif get(handles.kmedoidscluster, 'Value') == 1
    
    if size(data,1)==1
        labels=1;
    else
        k=min(k, size(data,1));
        labels = kmedoids(data', k)';
    end
end