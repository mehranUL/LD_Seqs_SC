function similarity = cosine_sim(a,b)
    if numel(a) ~= numel(b)
        error ( 'Size of two vectors are not the same!' );
    end
    similarity = dot(a,b)/(norm(a)*norm(b));
end