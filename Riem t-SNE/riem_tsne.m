function ydata = riem_tsne(X, no_dims, perplexity,M)
%Riem_tsne Performs Riemannian t-SNE [1,2,3] on dataset X on SPD Manifold 
% with the Log-Cholesky Riemannian metric [4]

% Input 

% Bib
% [1] - Bergsson, A., & Hauberg, S. (2022). Visualizing Riemannian 
% data with Rie-SNE. arXiv preprint arXiv:2203.09253.
% [2] - Laurens van der Maaten, 2010
% University of California, San Diego
% [3] - Van der Maaten, L., & Hinton, G. (2008). Visualizing data using
% t-SNE. Journal of machine learning research, 9(11).
% [4] - Lin, Z. (2019). Riemannian geometry of symmetric positive 
% definite matrices via Cholesky decomposition. SIAM Journal on 
% Matrix Analysis and Applications, 40(4), 1353-1370.



    if ~exist('labels', 'var')
        labels = [];
    end
    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
     if ~exist('initial_dims', 'var') || isempty(initial_dims)
        initial_dims = min(50, size(X, 2));
    end
    if ~exist('perplexity', 'var') || isempty(perplexity)
        perplexity = 30;
    end
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
        perplexity = initial_dims;
    else
        initial_solution = false;
    end
    
    
   % Compute pairwise distance matrix and volume ratios 
   [dim,~,N]=size(X);
   X=M.Manipulate.HPD2Chol(X);
   D=zeros(N,N);
   H=zeros(N,N);
   G=zeros(1,N);
   for i=1:N
       [~,g]=M.RiemOp.MetricTensor2(X(:,:,i));
       G(i)=g;
   end
   %--------------------------------
   for i=1:N
       for j=1:N
           if j<=i
               D(i,j)=M.RiemOp.Dist(X(:,:,i),X(:,:,j));
               H(i,j)=sqrt(G(i)/G(j));
           end
           H(i,j)=sqrt(G(i)/G(j));
       end        
       i
   end
   D=(D+D')/2;
   
    % Compute joint probabilities
    P = dh2p(D,H, perplexity, 1e-5);                                           % compute affinities using fixed perplexity
    clear D
    
    % Run t-SNE
    ydata = tsne_p(P,[],no_dims);

    



end



  
  