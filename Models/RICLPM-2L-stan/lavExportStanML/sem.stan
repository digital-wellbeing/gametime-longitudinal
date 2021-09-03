functions{

  vector[] sem_mean(vector[] alpha, real[,,] B, real[,,] gamma, int[] g, int k, int Ng, int gamind, real[,] meanx){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(gamind == 1){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (alpha[j] + to_matrix(gamma[,,j]) * to_vector(meanx[,j]));

      } else {
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * alpha[j];
      }
    }

    return evlv;
  }

  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;
    nidx = nov + nlvno0;

    iden = diag_matrix(rep_vector(1.0, k));

    if(nlvno0 > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:nidx){
          psivecinv[j,idx[i]] = 1/psi[idx[i],idx[i],j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,idx,idx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	  ldetcomp[j] = -2 * ldetcomp[j] + sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
	psimat[j] = to_matrix(psi[,,j]) + to_matrix(psi[,,j])' - diag_matrix(diagonal(to_matrix(psi[,,j])));

	ldetcomp[j] = log_determinant(psimat[j,idx,idx]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] - 2 * log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	}

	psimatinv[j] = psimat[j];
	psimatinv[j,1:nidx,1:nidx] = inverse_spd(psimat[j,idx,idx]);
        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,1:nidx,1:nidx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec[idx] - evlv[g[i],idx])' * siginv[g[i],1:nidx,1:nidx] * (xvec[idx] - evlv[g[i],idx]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }

  matrix fill_lower(matrix x){
    matrix[rows(x),cols(x)] newx;

    newx = x;
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        newx[j,i] = x[i,j];
      }
    }
    return newx;
  }
}

data{
  int N;
  int g[N];
  int lvind[8];
  int etaind[8];
  real sampmean[6,10];
  real meanx[6,10];
  int dummyov[2];
  int dummylv[2];
  vector[6] y[N];
  real lambdaframe[6,8,10];
  real thetaframe[6,6,10];
  real psiframe[8,8,10];
  real betaframe[8,8,10];
  real nuframe[6,1,10];
  real alphaframe[8,1,10];
  real lvrhoframe[8,8,10];
}

parameters{
  vector<lower=0>[80] psifree;
  vector[4] betafree;
  vector[60] nufree;
  vector<lower=0,upper=1>[40] lvrhofree;
  matrix[N, 8] etavec;
}

transformed parameters{
  real lambda[6,8,10];
  real theta[6,6,10];
  matrix[6,6] thetld[10];
  real psi[8,8,10];
  real beta[8,8,10];
  real nu[6,1,10];
  real alpha[8,1,10];
  real lvrho[8,8,10];
  real mu[N,6];
  matrix[N,8] eta;

  eta = rep_matrix(0, N, 8);

  lambda = lambdaframe;
  theta = thetaframe;
  psi = psiframe;
  beta = betaframe;
  nu = nuframe;
  alpha = alphaframe;
  lvrho = lvrhoframe;

  lambda[1,1,1] = 1;
  lambda[2,1,1] = 1;
  lambda[3,1,1] = 1;
  lambda[4,2,1] = 1;
  lambda[5,2,1] = 1;
  lambda[6,2,1] = 1;
  lambda[1,3,1] = 1;
  lambda[2,4,1] = 1;
  lambda[3,5,1] = 1;
  lambda[4,6,1] = 1;
  lambda[5,7,1] = 1;
  lambda[6,8,1] = 1;
  beta[4,3,1] = betafree[1];
  beta[4,6,1] = betafree[2];
  beta[7,3,1] = betafree[3];
  beta[7,6,1] = betafree[4];
  beta[5,4,1] = betafree[1];
  beta[5,7,1] = betafree[2];
  beta[8,4,1] = betafree[3];
  beta[8,7,1] = betafree[4];
  psi[1,1,1] = pow(psifree[1],2);
  psi[2,2,1] = pow(psifree[2],2);
  psi[3,3,1] = pow(psifree[3],2);
  psi[6,6,1] = pow(psifree[4],2);
  psi[4,4,1] = pow(psifree[5],2);
  psi[7,7,1] = pow(psifree[6],2);
  psi[5,5,1] = pow(psifree[7],2);
  psi[8,8,1] = pow(psifree[8],2);
  theta[1,1,1] = 0.001;
  theta[2,2,1] = 0.001;
  theta[3,3,1] = 0.001;
  theta[4,4,1] = 0.001;
  theta[5,5,1] = 0.001;
  theta[6,6,1] = 0.001;
  nu[1,1,1] = nufree[1];
  nu[2,1,1] = nufree[2];
  nu[3,1,1] = nufree[3];
  nu[4,1,1] = nufree[4];
  nu[5,1,1] = nufree[5];
  nu[6,1,1] = nufree[6];
  alpha[1,1,1] = 0;
  alpha[2,1,1] = 0;
  alpha[3,1,1] = 0;
  alpha[4,1,1] = 0;
  alpha[5,1,1] = 0;
  alpha[6,1,1] = 0;
  alpha[7,1,1] = 0;
  alpha[8,1,1] = 0;
  lvrho[3,6,1] = -1 + 2*lvrhofree[1];
  lvrho[4,7,1] = -1 + 2*lvrhofree[2];
  lvrho[5,8,1] = -1 + 2*lvrhofree[3];
  lvrho[1,2,1] = -1 + 2*lvrhofree[4];
  psi[3,6,1] = lvrho[3,6,1] * sqrt(psi[3,3,1] * psi[6,6,1]);
  psi[4,7,1] = lvrho[4,7,1] * sqrt(psi[4,4,1] * psi[7,7,1]);
  psi[5,8,1] = lvrho[5,8,1] * sqrt(psi[5,5,1] * psi[8,8,1]);
  psi[1,2,1] = lvrho[1,2,1] * sqrt(psi[1,1,1] * psi[2,2,1]);
  lambda[1,1,2] = 1;
  lambda[2,1,2] = 1;
  lambda[3,1,2] = 1;
  lambda[4,2,2] = 1;
  lambda[5,2,2] = 1;
  lambda[6,2,2] = 1;
  lambda[1,3,2] = 1;
  lambda[2,4,2] = 1;
  lambda[3,5,2] = 1;
  lambda[4,6,2] = 1;
  lambda[5,7,2] = 1;
  lambda[6,8,2] = 1;
  beta[4,3,2] = betafree[1];
  beta[4,6,2] = betafree[2];
  beta[7,3,2] = betafree[3];
  beta[7,6,2] = betafree[4];
  beta[5,4,2] = betafree[1];
  beta[5,7,2] = betafree[2];
  beta[8,4,2] = betafree[3];
  beta[8,7,2] = betafree[4];
  psi[1,1,2] = pow(psifree[9],2);
  psi[2,2,2] = pow(psifree[10],2);
  psi[3,3,2] = pow(psifree[11],2);
  psi[6,6,2] = pow(psifree[12],2);
  psi[4,4,2] = pow(psifree[13],2);
  psi[7,7,2] = pow(psifree[14],2);
  psi[5,5,2] = pow(psifree[15],2);
  psi[8,8,2] = pow(psifree[16],2);
  theta[1,1,2] = 0.001;
  theta[2,2,2] = 0.001;
  theta[3,3,2] = 0.001;
  theta[4,4,2] = 0.001;
  theta[5,5,2] = 0.001;
  theta[6,6,2] = 0.001;
  nu[1,1,2] = nufree[7];
  nu[2,1,2] = nufree[8];
  nu[3,1,2] = nufree[9];
  nu[4,1,2] = nufree[10];
  nu[5,1,2] = nufree[11];
  nu[6,1,2] = nufree[12];
  alpha[1,1,2] = 0;
  alpha[2,1,2] = 0;
  alpha[3,1,2] = 0;
  alpha[4,1,2] = 0;
  alpha[5,1,2] = 0;
  alpha[6,1,2] = 0;
  alpha[7,1,2] = 0;
  alpha[8,1,2] = 0;
  lvrho[3,6,2] = -1 + 2*lvrhofree[5];
  lvrho[4,7,2] = -1 + 2*lvrhofree[6];
  lvrho[5,8,2] = -1 + 2*lvrhofree[7];
  lvrho[1,2,2] = -1 + 2*lvrhofree[8];
  psi[3,6,2] = lvrho[3,6,2] * sqrt(psi[3,3,2] * psi[6,6,2]);
  psi[4,7,2] = lvrho[4,7,2] * sqrt(psi[4,4,2] * psi[7,7,2]);
  psi[5,8,2] = lvrho[5,8,2] * sqrt(psi[5,5,2] * psi[8,8,2]);
  psi[1,2,2] = lvrho[1,2,2] * sqrt(psi[1,1,2] * psi[2,2,2]);
  lambda[1,1,3] = 1;
  lambda[2,1,3] = 1;
  lambda[3,1,3] = 1;
  lambda[4,2,3] = 1;
  lambda[5,2,3] = 1;
  lambda[6,2,3] = 1;
  lambda[1,3,3] = 1;
  lambda[2,4,3] = 1;
  lambda[3,5,3] = 1;
  lambda[4,6,3] = 1;
  lambda[5,7,3] = 1;
  lambda[6,8,3] = 1;
  beta[4,3,3] = betafree[1];
  beta[4,6,3] = betafree[2];
  beta[7,3,3] = betafree[3];
  beta[7,6,3] = betafree[4];
  beta[5,4,3] = betafree[1];
  beta[5,7,3] = betafree[2];
  beta[8,4,3] = betafree[3];
  beta[8,7,3] = betafree[4];
  psi[1,1,3] = pow(psifree[17],2);
  psi[2,2,3] = pow(psifree[18],2);
  psi[3,3,3] = pow(psifree[19],2);
  psi[6,6,3] = pow(psifree[20],2);
  psi[4,4,3] = pow(psifree[21],2);
  psi[7,7,3] = pow(psifree[22],2);
  psi[5,5,3] = pow(psifree[23],2);
  psi[8,8,3] = pow(psifree[24],2);
  theta[1,1,3] = 0.001;
  theta[2,2,3] = 0.001;
  theta[3,3,3] = 0.001;
  theta[4,4,3] = 0.001;
  theta[5,5,3] = 0.001;
  theta[6,6,3] = 0.001;
  nu[1,1,3] = nufree[13];
  nu[2,1,3] = nufree[14];
  nu[3,1,3] = nufree[15];
  nu[4,1,3] = nufree[16];
  nu[5,1,3] = nufree[17];
  nu[6,1,3] = nufree[18];
  alpha[1,1,3] = 0;
  alpha[2,1,3] = 0;
  alpha[3,1,3] = 0;
  alpha[4,1,3] = 0;
  alpha[5,1,3] = 0;
  alpha[6,1,3] = 0;
  alpha[7,1,3] = 0;
  alpha[8,1,3] = 0;
  lvrho[3,6,3] = -1 + 2*lvrhofree[9];
  lvrho[4,7,3] = -1 + 2*lvrhofree[10];
  lvrho[5,8,3] = -1 + 2*lvrhofree[11];
  lvrho[1,2,3] = -1 + 2*lvrhofree[12];
  psi[3,6,3] = lvrho[3,6,3] * sqrt(psi[3,3,3] * psi[6,6,3]);
  psi[4,7,3] = lvrho[4,7,3] * sqrt(psi[4,4,3] * psi[7,7,3]);
  psi[5,8,3] = lvrho[5,8,3] * sqrt(psi[5,5,3] * psi[8,8,3]);
  psi[1,2,3] = lvrho[1,2,3] * sqrt(psi[1,1,3] * psi[2,2,3]);
  lambda[1,1,4] = 1;
  lambda[2,1,4] = 1;
  lambda[3,1,4] = 1;
  lambda[4,2,4] = 1;
  lambda[5,2,4] = 1;
  lambda[6,2,4] = 1;
  lambda[1,3,4] = 1;
  lambda[2,4,4] = 1;
  lambda[3,5,4] = 1;
  lambda[4,6,4] = 1;
  lambda[5,7,4] = 1;
  lambda[6,8,4] = 1;
  beta[4,3,4] = betafree[1];
  beta[4,6,4] = betafree[2];
  beta[7,3,4] = betafree[3];
  beta[7,6,4] = betafree[4];
  beta[5,4,4] = betafree[1];
  beta[5,7,4] = betafree[2];
  beta[8,4,4] = betafree[3];
  beta[8,7,4] = betafree[4];
  psi[1,1,4] = pow(psifree[25],2);
  psi[2,2,4] = pow(psifree[26],2);
  psi[3,3,4] = pow(psifree[27],2);
  psi[6,6,4] = pow(psifree[28],2);
  psi[4,4,4] = pow(psifree[29],2);
  psi[7,7,4] = pow(psifree[30],2);
  psi[5,5,4] = pow(psifree[31],2);
  psi[8,8,4] = pow(psifree[32],2);
  theta[1,1,4] = 0.001;
  theta[2,2,4] = 0.001;
  theta[3,3,4] = 0.001;
  theta[4,4,4] = 0.001;
  theta[5,5,4] = 0.001;
  theta[6,6,4] = 0.001;
  nu[1,1,4] = nufree[19];
  nu[2,1,4] = nufree[20];
  nu[3,1,4] = nufree[21];
  nu[4,1,4] = nufree[22];
  nu[5,1,4] = nufree[23];
  nu[6,1,4] = nufree[24];
  alpha[1,1,4] = 0;
  alpha[2,1,4] = 0;
  alpha[3,1,4] = 0;
  alpha[4,1,4] = 0;
  alpha[5,1,4] = 0;
  alpha[6,1,4] = 0;
  alpha[7,1,4] = 0;
  alpha[8,1,4] = 0;
  lvrho[3,6,4] = -1 + 2*lvrhofree[13];
  lvrho[4,7,4] = -1 + 2*lvrhofree[14];
  lvrho[5,8,4] = -1 + 2*lvrhofree[15];
  lvrho[1,2,4] = -1 + 2*lvrhofree[16];
  psi[3,6,4] = lvrho[3,6,4] * sqrt(psi[3,3,4] * psi[6,6,4]);
  psi[4,7,4] = lvrho[4,7,4] * sqrt(psi[4,4,4] * psi[7,7,4]);
  psi[5,8,4] = lvrho[5,8,4] * sqrt(psi[5,5,4] * psi[8,8,4]);
  psi[1,2,4] = lvrho[1,2,4] * sqrt(psi[1,1,4] * psi[2,2,4]);
  lambda[1,1,5] = 1;
  lambda[2,1,5] = 1;
  lambda[3,1,5] = 1;
  lambda[4,2,5] = 1;
  lambda[5,2,5] = 1;
  lambda[6,2,5] = 1;
  lambda[1,3,5] = 1;
  lambda[2,4,5] = 1;
  lambda[3,5,5] = 1;
  lambda[4,6,5] = 1;
  lambda[5,7,5] = 1;
  lambda[6,8,5] = 1;
  beta[4,3,5] = betafree[1];
  beta[4,6,5] = betafree[2];
  beta[7,3,5] = betafree[3];
  beta[7,6,5] = betafree[4];
  beta[5,4,5] = betafree[1];
  beta[5,7,5] = betafree[2];
  beta[8,4,5] = betafree[3];
  beta[8,7,5] = betafree[4];
  psi[1,1,5] = pow(psifree[33],2);
  psi[2,2,5] = pow(psifree[34],2);
  psi[3,3,5] = pow(psifree[35],2);
  psi[6,6,5] = pow(psifree[36],2);
  psi[4,4,5] = pow(psifree[37],2);
  psi[7,7,5] = pow(psifree[38],2);
  psi[5,5,5] = pow(psifree[39],2);
  psi[8,8,5] = pow(psifree[40],2);
  theta[1,1,5] = 0.001;
  theta[2,2,5] = 0.001;
  theta[3,3,5] = 0.001;
  theta[4,4,5] = 0.001;
  theta[5,5,5] = 0.001;
  theta[6,6,5] = 0.001;
  nu[1,1,5] = nufree[25];
  nu[2,1,5] = nufree[26];
  nu[3,1,5] = nufree[27];
  nu[4,1,5] = nufree[28];
  nu[5,1,5] = nufree[29];
  nu[6,1,5] = nufree[30];
  alpha[1,1,5] = 0;
  alpha[2,1,5] = 0;
  alpha[3,1,5] = 0;
  alpha[4,1,5] = 0;
  alpha[5,1,5] = 0;
  alpha[6,1,5] = 0;
  alpha[7,1,5] = 0;
  alpha[8,1,5] = 0;
  lvrho[3,6,5] = -1 + 2*lvrhofree[17];
  lvrho[4,7,5] = -1 + 2*lvrhofree[18];
  lvrho[5,8,5] = -1 + 2*lvrhofree[19];
  lvrho[1,2,5] = -1 + 2*lvrhofree[20];
  psi[3,6,5] = lvrho[3,6,5] * sqrt(psi[3,3,5] * psi[6,6,5]);
  psi[4,7,5] = lvrho[4,7,5] * sqrt(psi[4,4,5] * psi[7,7,5]);
  psi[5,8,5] = lvrho[5,8,5] * sqrt(psi[5,5,5] * psi[8,8,5]);
  psi[1,2,5] = lvrho[1,2,5] * sqrt(psi[1,1,5] * psi[2,2,5]);
  lambda[1,1,6] = 1;
  lambda[2,1,6] = 1;
  lambda[3,1,6] = 1;
  lambda[4,2,6] = 1;
  lambda[5,2,6] = 1;
  lambda[6,2,6] = 1;
  lambda[1,3,6] = 1;
  lambda[2,4,6] = 1;
  lambda[3,5,6] = 1;
  lambda[4,6,6] = 1;
  lambda[5,7,6] = 1;
  lambda[6,8,6] = 1;
  beta[4,3,6] = betafree[1];
  beta[4,6,6] = betafree[2];
  beta[7,3,6] = betafree[3];
  beta[7,6,6] = betafree[4];
  beta[5,4,6] = betafree[1];
  beta[5,7,6] = betafree[2];
  beta[8,4,6] = betafree[3];
  beta[8,7,6] = betafree[4];
  psi[1,1,6] = pow(psifree[41],2);
  psi[2,2,6] = pow(psifree[42],2);
  psi[3,3,6] = pow(psifree[43],2);
  psi[6,6,6] = pow(psifree[44],2);
  psi[4,4,6] = pow(psifree[45],2);
  psi[7,7,6] = pow(psifree[46],2);
  psi[5,5,6] = pow(psifree[47],2);
  psi[8,8,6] = pow(psifree[48],2);
  theta[1,1,6] = 0.001;
  theta[2,2,6] = 0.001;
  theta[3,3,6] = 0.001;
  theta[4,4,6] = 0.001;
  theta[5,5,6] = 0.001;
  theta[6,6,6] = 0.001;
  nu[1,1,6] = nufree[31];
  nu[2,1,6] = nufree[32];
  nu[3,1,6] = nufree[33];
  nu[4,1,6] = nufree[34];
  nu[5,1,6] = nufree[35];
  nu[6,1,6] = nufree[36];
  alpha[1,1,6] = 0;
  alpha[2,1,6] = 0;
  alpha[3,1,6] = 0;
  alpha[4,1,6] = 0;
  alpha[5,1,6] = 0;
  alpha[6,1,6] = 0;
  alpha[7,1,6] = 0;
  alpha[8,1,6] = 0;
  lvrho[3,6,6] = -1 + 2*lvrhofree[21];
  lvrho[4,7,6] = -1 + 2*lvrhofree[22];
  lvrho[5,8,6] = -1 + 2*lvrhofree[23];
  lvrho[1,2,6] = -1 + 2*lvrhofree[24];
  psi[3,6,6] = lvrho[3,6,6] * sqrt(psi[3,3,6] * psi[6,6,6]);
  psi[4,7,6] = lvrho[4,7,6] * sqrt(psi[4,4,6] * psi[7,7,6]);
  psi[5,8,6] = lvrho[5,8,6] * sqrt(psi[5,5,6] * psi[8,8,6]);
  psi[1,2,6] = lvrho[1,2,6] * sqrt(psi[1,1,6] * psi[2,2,6]);
  lambda[1,1,7] = 1;
  lambda[2,1,7] = 1;
  lambda[3,1,7] = 1;
  lambda[4,2,7] = 1;
  lambda[5,2,7] = 1;
  lambda[6,2,7] = 1;
  lambda[1,3,7] = 1;
  lambda[2,4,7] = 1;
  lambda[3,5,7] = 1;
  lambda[4,6,7] = 1;
  lambda[5,7,7] = 1;
  lambda[6,8,7] = 1;
  beta[4,3,7] = betafree[1];
  beta[4,6,7] = betafree[2];
  beta[7,3,7] = betafree[3];
  beta[7,6,7] = betafree[4];
  beta[5,4,7] = betafree[1];
  beta[5,7,7] = betafree[2];
  beta[8,4,7] = betafree[3];
  beta[8,7,7] = betafree[4];
  psi[1,1,7] = pow(psifree[49],2);
  psi[2,2,7] = pow(psifree[50],2);
  psi[3,3,7] = pow(psifree[51],2);
  psi[6,6,7] = pow(psifree[52],2);
  psi[4,4,7] = pow(psifree[53],2);
  psi[7,7,7] = pow(psifree[54],2);
  psi[5,5,7] = pow(psifree[55],2);
  psi[8,8,7] = pow(psifree[56],2);
  theta[1,1,7] = 0.001;
  theta[2,2,7] = 0.001;
  theta[3,3,7] = 0.001;
  theta[4,4,7] = 0.001;
  theta[5,5,7] = 0.001;
  theta[6,6,7] = 0.001;
  nu[1,1,7] = nufree[37];
  nu[2,1,7] = nufree[38];
  nu[3,1,7] = nufree[39];
  nu[4,1,7] = nufree[40];
  nu[5,1,7] = nufree[41];
  nu[6,1,7] = nufree[42];
  alpha[1,1,7] = 0;
  alpha[2,1,7] = 0;
  alpha[3,1,7] = 0;
  alpha[4,1,7] = 0;
  alpha[5,1,7] = 0;
  alpha[6,1,7] = 0;
  alpha[7,1,7] = 0;
  alpha[8,1,7] = 0;
  lvrho[3,6,7] = -1 + 2*lvrhofree[25];
  lvrho[4,7,7] = -1 + 2*lvrhofree[26];
  lvrho[5,8,7] = -1 + 2*lvrhofree[27];
  lvrho[1,2,7] = -1 + 2*lvrhofree[28];
  psi[3,6,7] = lvrho[3,6,7] * sqrt(psi[3,3,7] * psi[6,6,7]);
  psi[4,7,7] = lvrho[4,7,7] * sqrt(psi[4,4,7] * psi[7,7,7]);
  psi[5,8,7] = lvrho[5,8,7] * sqrt(psi[5,5,7] * psi[8,8,7]);
  psi[1,2,7] = lvrho[1,2,7] * sqrt(psi[1,1,7] * psi[2,2,7]);
  lambda[1,1,8] = 1;
  lambda[2,1,8] = 1;
  lambda[3,1,8] = 1;
  lambda[4,2,8] = 1;
  lambda[5,2,8] = 1;
  lambda[6,2,8] = 1;
  lambda[1,3,8] = 1;
  lambda[2,4,8] = 1;
  lambda[3,5,8] = 1;
  lambda[4,6,8] = 1;
  lambda[5,7,8] = 1;
  lambda[6,8,8] = 1;
  beta[4,3,8] = betafree[1];
  beta[4,6,8] = betafree[2];
  beta[7,3,8] = betafree[3];
  beta[7,6,8] = betafree[4];
  beta[5,4,8] = betafree[1];
  beta[5,7,8] = betafree[2];
  beta[8,4,8] = betafree[3];
  beta[8,7,8] = betafree[4];
  psi[1,1,8] = pow(psifree[57],2);
  psi[2,2,8] = pow(psifree[58],2);
  psi[3,3,8] = pow(psifree[59],2);
  psi[6,6,8] = pow(psifree[60],2);
  psi[4,4,8] = pow(psifree[61],2);
  psi[7,7,8] = pow(psifree[62],2);
  psi[5,5,8] = pow(psifree[63],2);
  psi[8,8,8] = pow(psifree[64],2);
  theta[1,1,8] = 0.001;
  theta[2,2,8] = 0.001;
  theta[3,3,8] = 0.001;
  theta[4,4,8] = 0.001;
  theta[5,5,8] = 0.001;
  theta[6,6,8] = 0.001;
  nu[1,1,8] = nufree[43];
  nu[2,1,8] = nufree[44];
  nu[3,1,8] = nufree[45];
  nu[4,1,8] = nufree[46];
  nu[5,1,8] = nufree[47];
  nu[6,1,8] = nufree[48];
  alpha[1,1,8] = 0;
  alpha[2,1,8] = 0;
  alpha[3,1,8] = 0;
  alpha[4,1,8] = 0;
  alpha[5,1,8] = 0;
  alpha[6,1,8] = 0;
  alpha[7,1,8] = 0;
  alpha[8,1,8] = 0;
  lvrho[3,6,8] = -1 + 2*lvrhofree[29];
  lvrho[4,7,8] = -1 + 2*lvrhofree[30];
  lvrho[5,8,8] = -1 + 2*lvrhofree[31];
  lvrho[1,2,8] = -1 + 2*lvrhofree[32];
  psi[3,6,8] = lvrho[3,6,8] * sqrt(psi[3,3,8] * psi[6,6,8]);
  psi[4,7,8] = lvrho[4,7,8] * sqrt(psi[4,4,8] * psi[7,7,8]);
  psi[5,8,8] = lvrho[5,8,8] * sqrt(psi[5,5,8] * psi[8,8,8]);
  psi[1,2,8] = lvrho[1,2,8] * sqrt(psi[1,1,8] * psi[2,2,8]);
  lambda[1,1,9] = 1;
  lambda[2,1,9] = 1;
  lambda[3,1,9] = 1;
  lambda[4,2,9] = 1;
  lambda[5,2,9] = 1;
  lambda[6,2,9] = 1;
  lambda[1,3,9] = 1;
  lambda[2,4,9] = 1;
  lambda[3,5,9] = 1;
  lambda[4,6,9] = 1;
  lambda[5,7,9] = 1;
  lambda[6,8,9] = 1;
  beta[4,3,9] = betafree[1];
  beta[4,6,9] = betafree[2];
  beta[7,3,9] = betafree[3];
  beta[7,6,9] = betafree[4];
  beta[5,4,9] = betafree[1];
  beta[5,7,9] = betafree[2];
  beta[8,4,9] = betafree[3];
  beta[8,7,9] = betafree[4];
  psi[1,1,9] = pow(psifree[65],2);
  psi[2,2,9] = pow(psifree[66],2);
  psi[3,3,9] = pow(psifree[67],2);
  psi[6,6,9] = pow(psifree[68],2);
  psi[4,4,9] = pow(psifree[69],2);
  psi[7,7,9] = pow(psifree[70],2);
  psi[5,5,9] = pow(psifree[71],2);
  psi[8,8,9] = pow(psifree[72],2);
  theta[1,1,9] = 0.001;
  theta[2,2,9] = 0.001;
  theta[3,3,9] = 0.001;
  theta[4,4,9] = 0.001;
  theta[5,5,9] = 0.001;
  theta[6,6,9] = 0.001;
  nu[1,1,9] = nufree[49];
  nu[2,1,9] = nufree[50];
  nu[3,1,9] = nufree[51];
  nu[4,1,9] = nufree[52];
  nu[5,1,9] = nufree[53];
  nu[6,1,9] = nufree[54];
  alpha[1,1,9] = 0;
  alpha[2,1,9] = 0;
  alpha[3,1,9] = 0;
  alpha[4,1,9] = 0;
  alpha[5,1,9] = 0;
  alpha[6,1,9] = 0;
  alpha[7,1,9] = 0;
  alpha[8,1,9] = 0;
  lvrho[3,6,9] = -1 + 2*lvrhofree[33];
  lvrho[4,7,9] = -1 + 2*lvrhofree[34];
  lvrho[5,8,9] = -1 + 2*lvrhofree[35];
  lvrho[1,2,9] = -1 + 2*lvrhofree[36];
  psi[3,6,9] = lvrho[3,6,9] * sqrt(psi[3,3,9] * psi[6,6,9]);
  psi[4,7,9] = lvrho[4,7,9] * sqrt(psi[4,4,9] * psi[7,7,9]);
  psi[5,8,9] = lvrho[5,8,9] * sqrt(psi[5,5,9] * psi[8,8,9]);
  psi[1,2,9] = lvrho[1,2,9] * sqrt(psi[1,1,9] * psi[2,2,9]);
  lambda[1,1,10] = 1;
  lambda[2,1,10] = 1;
  lambda[3,1,10] = 1;
  lambda[4,2,10] = 1;
  lambda[5,2,10] = 1;
  lambda[6,2,10] = 1;
  lambda[1,3,10] = 1;
  lambda[2,4,10] = 1;
  lambda[3,5,10] = 1;
  lambda[4,6,10] = 1;
  lambda[5,7,10] = 1;
  lambda[6,8,10] = 1;
  beta[4,3,10] = betafree[1];
  beta[4,6,10] = betafree[2];
  beta[7,3,10] = betafree[3];
  beta[7,6,10] = betafree[4];
  beta[5,4,10] = betafree[1];
  beta[5,7,10] = betafree[2];
  beta[8,4,10] = betafree[3];
  beta[8,7,10] = betafree[4];
  psi[1,1,10] = pow(psifree[73],2);
  psi[2,2,10] = pow(psifree[74],2);
  psi[3,3,10] = pow(psifree[75],2);
  psi[6,6,10] = pow(psifree[76],2);
  psi[4,4,10] = pow(psifree[77],2);
  psi[7,7,10] = pow(psifree[78],2);
  psi[5,5,10] = pow(psifree[79],2);
  psi[8,8,10] = pow(psifree[80],2);
  theta[1,1,10] = 0.001;
  theta[2,2,10] = 0.001;
  theta[3,3,10] = 0.001;
  theta[4,4,10] = 0.001;
  theta[5,5,10] = 0.001;
  theta[6,6,10] = 0.001;
  nu[1,1,10] = nufree[55];
  nu[2,1,10] = nufree[56];
  nu[3,1,10] = nufree[57];
  nu[4,1,10] = nufree[58];
  nu[5,1,10] = nufree[59];
  nu[6,1,10] = nufree[60];
  alpha[1,1,10] = 0;
  alpha[2,1,10] = 0;
  alpha[3,1,10] = 0;
  alpha[4,1,10] = 0;
  alpha[5,1,10] = 0;
  alpha[6,1,10] = 0;
  alpha[7,1,10] = 0;
  alpha[8,1,10] = 0;
  lvrho[3,6,10] = -1 + 2*lvrhofree[37];
  lvrho[4,7,10] = -1 + 2*lvrhofree[38];
  lvrho[5,8,10] = -1 + 2*lvrhofree[39];
  lvrho[1,2,10] = -1 + 2*lvrhofree[40];
  psi[3,6,10] = lvrho[3,6,10] * sqrt(psi[3,3,10] * psi[6,6,10]);
  psi[4,7,10] = lvrho[4,7,10] * sqrt(psi[4,4,10] * psi[7,7,10]);
  psi[5,8,10] = lvrho[5,8,10] * sqrt(psi[5,5,10] * psi[8,8,10]);
  psi[1,2,10] = lvrho[1,2,10] * sqrt(psi[1,1,10] * psi[2,2,10]);

  // mu definitions
  for(i in 1:N) {
    eta[i,1:8] = etavec[i];

    mu[i,1] = nu[1,1,g[i]] + lambda[1,1,g[i]]*eta[i,1] + lambda[1,3,g[i]]*eta[i,3];
    mu[i,2] = nu[2,1,g[i]] + lambda[2,1,g[i]]*eta[i,1] + lambda[2,4,g[i]]*eta[i,4];
    mu[i,3] = nu[3,1,g[i]] + lambda[3,1,g[i]]*eta[i,1] + lambda[3,5,g[i]]*eta[i,5];
    mu[i,4] = nu[4,1,g[i]] + lambda[4,2,g[i]]*eta[i,2] + lambda[4,6,g[i]]*eta[i,6];
    mu[i,5] = nu[5,1,g[i]] + lambda[5,2,g[i]]*eta[i,2] + lambda[5,7,g[i]]*eta[i,7];
    mu[i,6] = nu[6,1,g[i]] + lambda[6,2,g[i]]*eta[i,2] + lambda[6,8,g[i]]*eta[i,8];
  }

  for(j in 1:10){
    thetld[j] = fill_lower(to_matrix(theta[,,j]));
    thetld[j] = cholesky_decompose(thetld[j]);
  }

}

model {
  for(i in 1:N) {
    y[i] ~ multi_normal_cholesky(to_vector(mu[i,1:6]), thetld[g[i]]);
  }

  eta ~ sem_lv(alpha, beta, psi, beta, 0, meanx, g, 8, N, 10, 0, 1, 8, etaind, 8);

  // Priors
  target += normal_lpdf(betafree[1] | 0,10);
  target += normal_lpdf(betafree[2] | 0,10);
  target += normal_lpdf(betafree[3] | 0,10);
  target += normal_lpdf(betafree[4] | 0,10);
  target += gamma_lpdf(psifree[1] | 1,.5);
  target += gamma_lpdf(psifree[2] | 1,.5);
  target += gamma_lpdf(psifree[3] | 1,.5);
  target += gamma_lpdf(psifree[4] | 1,.5);
  target += gamma_lpdf(psifree[5] | 1,.5);
  target += gamma_lpdf(psifree[6] | 1,.5);
  target += gamma_lpdf(psifree[7] | 1,.5);
  target += gamma_lpdf(psifree[8] | 1,.5);
  target += normal_lpdf(nufree[1] | 0,32);
  target += normal_lpdf(nufree[2] | 0,32);
  target += normal_lpdf(nufree[3] | 0,32);
  target += normal_lpdf(nufree[4] | 0,32);
  target += normal_lpdf(nufree[5] | 0,32);
  target += normal_lpdf(nufree[6] | 0,32);
  target += beta_lpdf(lvrhofree[1] | 1,1);
  target += beta_lpdf(lvrhofree[2] | 1,1);
  target += beta_lpdf(lvrhofree[3] | 1,1);
  target += beta_lpdf(lvrhofree[4] | 1,1);
  target += gamma_lpdf(psifree[9] | 1,.5);
  target += gamma_lpdf(psifree[10] | 1,.5);
  target += gamma_lpdf(psifree[11] | 1,.5);
  target += gamma_lpdf(psifree[12] | 1,.5);
  target += gamma_lpdf(psifree[13] | 1,.5);
  target += gamma_lpdf(psifree[14] | 1,.5);
  target += gamma_lpdf(psifree[15] | 1,.5);
  target += gamma_lpdf(psifree[16] | 1,.5);
  target += normal_lpdf(nufree[7] | 0,32);
  target += normal_lpdf(nufree[8] | 0,32);
  target += normal_lpdf(nufree[9] | 0,32);
  target += normal_lpdf(nufree[10] | 0,32);
  target += normal_lpdf(nufree[11] | 0,32);
  target += normal_lpdf(nufree[12] | 0,32);
  target += beta_lpdf(lvrhofree[5] | 1,1);
  target += beta_lpdf(lvrhofree[6] | 1,1);
  target += beta_lpdf(lvrhofree[7] | 1,1);
  target += beta_lpdf(lvrhofree[8] | 1,1);
  target += gamma_lpdf(psifree[17] | 1,.5);
  target += gamma_lpdf(psifree[18] | 1,.5);
  target += gamma_lpdf(psifree[19] | 1,.5);
  target += gamma_lpdf(psifree[20] | 1,.5);
  target += gamma_lpdf(psifree[21] | 1,.5);
  target += gamma_lpdf(psifree[22] | 1,.5);
  target += gamma_lpdf(psifree[23] | 1,.5);
  target += gamma_lpdf(psifree[24] | 1,.5);
  target += normal_lpdf(nufree[13] | 0,32);
  target += normal_lpdf(nufree[14] | 0,32);
  target += normal_lpdf(nufree[15] | 0,32);
  target += normal_lpdf(nufree[16] | 0,32);
  target += normal_lpdf(nufree[17] | 0,32);
  target += normal_lpdf(nufree[18] | 0,32);
  target += beta_lpdf(lvrhofree[9] | 1,1);
  target += beta_lpdf(lvrhofree[10] | 1,1);
  target += beta_lpdf(lvrhofree[11] | 1,1);
  target += beta_lpdf(lvrhofree[12] | 1,1);
  target += gamma_lpdf(psifree[25] | 1,.5);
  target += gamma_lpdf(psifree[26] | 1,.5);
  target += gamma_lpdf(psifree[27] | 1,.5);
  target += gamma_lpdf(psifree[28] | 1,.5);
  target += gamma_lpdf(psifree[29] | 1,.5);
  target += gamma_lpdf(psifree[30] | 1,.5);
  target += gamma_lpdf(psifree[31] | 1,.5);
  target += gamma_lpdf(psifree[32] | 1,.5);
  target += normal_lpdf(nufree[19] | 0,32);
  target += normal_lpdf(nufree[20] | 0,32);
  target += normal_lpdf(nufree[21] | 0,32);
  target += normal_lpdf(nufree[22] | 0,32);
  target += normal_lpdf(nufree[23] | 0,32);
  target += normal_lpdf(nufree[24] | 0,32);
  target += beta_lpdf(lvrhofree[13] | 1,1);
  target += beta_lpdf(lvrhofree[14] | 1,1);
  target += beta_lpdf(lvrhofree[15] | 1,1);
  target += beta_lpdf(lvrhofree[16] | 1,1);
  target += gamma_lpdf(psifree[33] | 1,.5);
  target += gamma_lpdf(psifree[34] | 1,.5);
  target += gamma_lpdf(psifree[35] | 1,.5);
  target += gamma_lpdf(psifree[36] | 1,.5);
  target += gamma_lpdf(psifree[37] | 1,.5);
  target += gamma_lpdf(psifree[38] | 1,.5);
  target += gamma_lpdf(psifree[39] | 1,.5);
  target += gamma_lpdf(psifree[40] | 1,.5);
  target += normal_lpdf(nufree[25] | 0,32);
  target += normal_lpdf(nufree[26] | 0,32);
  target += normal_lpdf(nufree[27] | 0,32);
  target += normal_lpdf(nufree[28] | 0,32);
  target += normal_lpdf(nufree[29] | 0,32);
  target += normal_lpdf(nufree[30] | 0,32);
  target += beta_lpdf(lvrhofree[17] | 1,1);
  target += beta_lpdf(lvrhofree[18] | 1,1);
  target += beta_lpdf(lvrhofree[19] | 1,1);
  target += beta_lpdf(lvrhofree[20] | 1,1);
  target += gamma_lpdf(psifree[41] | 1,.5);
  target += gamma_lpdf(psifree[42] | 1,.5);
  target += gamma_lpdf(psifree[43] | 1,.5);
  target += gamma_lpdf(psifree[44] | 1,.5);
  target += gamma_lpdf(psifree[45] | 1,.5);
  target += gamma_lpdf(psifree[46] | 1,.5);
  target += gamma_lpdf(psifree[47] | 1,.5);
  target += gamma_lpdf(psifree[48] | 1,.5);
  target += normal_lpdf(nufree[31] | 0,32);
  target += normal_lpdf(nufree[32] | 0,32);
  target += normal_lpdf(nufree[33] | 0,32);
  target += normal_lpdf(nufree[34] | 0,32);
  target += normal_lpdf(nufree[35] | 0,32);
  target += normal_lpdf(nufree[36] | 0,32);
  target += beta_lpdf(lvrhofree[21] | 1,1);
  target += beta_lpdf(lvrhofree[22] | 1,1);
  target += beta_lpdf(lvrhofree[23] | 1,1);
  target += beta_lpdf(lvrhofree[24] | 1,1);
  target += gamma_lpdf(psifree[49] | 1,.5);
  target += gamma_lpdf(psifree[50] | 1,.5);
  target += gamma_lpdf(psifree[51] | 1,.5);
  target += gamma_lpdf(psifree[52] | 1,.5);
  target += gamma_lpdf(psifree[53] | 1,.5);
  target += gamma_lpdf(psifree[54] | 1,.5);
  target += gamma_lpdf(psifree[55] | 1,.5);
  target += gamma_lpdf(psifree[56] | 1,.5);
  target += normal_lpdf(nufree[37] | 0,32);
  target += normal_lpdf(nufree[38] | 0,32);
  target += normal_lpdf(nufree[39] | 0,32);
  target += normal_lpdf(nufree[40] | 0,32);
  target += normal_lpdf(nufree[41] | 0,32);
  target += normal_lpdf(nufree[42] | 0,32);
  target += beta_lpdf(lvrhofree[25] | 1,1);
  target += beta_lpdf(lvrhofree[26] | 1,1);
  target += beta_lpdf(lvrhofree[27] | 1,1);
  target += beta_lpdf(lvrhofree[28] | 1,1);
  target += gamma_lpdf(psifree[57] | 1,.5);
  target += gamma_lpdf(psifree[58] | 1,.5);
  target += gamma_lpdf(psifree[59] | 1,.5);
  target += gamma_lpdf(psifree[60] | 1,.5);
  target += gamma_lpdf(psifree[61] | 1,.5);
  target += gamma_lpdf(psifree[62] | 1,.5);
  target += gamma_lpdf(psifree[63] | 1,.5);
  target += gamma_lpdf(psifree[64] | 1,.5);
  target += normal_lpdf(nufree[43] | 0,32);
  target += normal_lpdf(nufree[44] | 0,32);
  target += normal_lpdf(nufree[45] | 0,32);
  target += normal_lpdf(nufree[46] | 0,32);
  target += normal_lpdf(nufree[47] | 0,32);
  target += normal_lpdf(nufree[48] | 0,32);
  target += beta_lpdf(lvrhofree[29] | 1,1);
  target += beta_lpdf(lvrhofree[30] | 1,1);
  target += beta_lpdf(lvrhofree[31] | 1,1);
  target += beta_lpdf(lvrhofree[32] | 1,1);
  target += gamma_lpdf(psifree[65] | 1,.5);
  target += gamma_lpdf(psifree[66] | 1,.5);
  target += gamma_lpdf(psifree[67] | 1,.5);
  target += gamma_lpdf(psifree[68] | 1,.5);
  target += gamma_lpdf(psifree[69] | 1,.5);
  target += gamma_lpdf(psifree[70] | 1,.5);
  target += gamma_lpdf(psifree[71] | 1,.5);
  target += gamma_lpdf(psifree[72] | 1,.5);
  target += normal_lpdf(nufree[49] | 0,32);
  target += normal_lpdf(nufree[50] | 0,32);
  target += normal_lpdf(nufree[51] | 0,32);
  target += normal_lpdf(nufree[52] | 0,32);
  target += normal_lpdf(nufree[53] | 0,32);
  target += normal_lpdf(nufree[54] | 0,32);
  target += beta_lpdf(lvrhofree[33] | 1,1);
  target += beta_lpdf(lvrhofree[34] | 1,1);
  target += beta_lpdf(lvrhofree[35] | 1,1);
  target += beta_lpdf(lvrhofree[36] | 1,1);
  target += gamma_lpdf(psifree[73] | 1,.5);
  target += gamma_lpdf(psifree[74] | 1,.5);
  target += gamma_lpdf(psifree[75] | 1,.5);
  target += gamma_lpdf(psifree[76] | 1,.5);
  target += gamma_lpdf(psifree[77] | 1,.5);
  target += gamma_lpdf(psifree[78] | 1,.5);
  target += gamma_lpdf(psifree[79] | 1,.5);
  target += gamma_lpdf(psifree[80] | 1,.5);
  target += normal_lpdf(nufree[55] | 0,32);
  target += normal_lpdf(nufree[56] | 0,32);
  target += normal_lpdf(nufree[57] | 0,32);
  target += normal_lpdf(nufree[58] | 0,32);
  target += normal_lpdf(nufree[59] | 0,32);
  target += normal_lpdf(nufree[60] | 0,32);
  target += beta_lpdf(lvrhofree[37] | 1,1);
  target += beta_lpdf(lvrhofree[38] | 1,1);
  target += beta_lpdf(lvrhofree[39] | 1,1);
  target += beta_lpdf(lvrhofree[40] | 1,1);
}