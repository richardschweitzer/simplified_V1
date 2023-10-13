get_log_gabor_heiko <- function(imSize, # size of image in pixels
                                degSize, # size of image in degrees of visual angle
                                freq,
                                orientation,
                                bw=c(0.5945, 0.2965),
                                phase,
                                scaler=1) {
  # translated to R from Heiko Schuett's Matlab function:
  # https://github.com/wichmann-lab/spatial-vision-model/blob/master/image_logGabor.m
  # Richard made some additional comments
  
  # Heiko says:
  # function image = image_logGabor(imSize,degSize,freq,contrast,center,orientation,L,bw,phase)
  # this function creates a Gabor patch with given orientation, frequency and
  # size (rad and deg of visual angle) in the center of the image
  #
  # 0 phase is the odd positive phase
  #
  # numerical normalization to 1 at peak of the Gabor
  
  require(pracma)
  
  # scale?
  if (scaler!=1) {
    imSize <- imSize * scaler
  }
  
  # make image, get frequencies
  scr_ppd <- imSize / degSize
  x = 1:imSize[2]
  y = 1:imSize[1]
  x_spat <- (x-ceil(mean(x))) / scr_ppd[2] # for later
  y_spat <- (y-ceil(mean(y))) / scr_ppd[1]
  x = (x-ceil(mean(x))) / degSize[2]
  y = (y-ceil(mean(y))) / degSize[1]
  
  # make frequency and orientation space
  meshi = meshgrid(x, y)
  f = sqrt(meshi$X^2+meshi$Y^2)
  orient0 = atan2(meshi$Y, meshi$X)
  
  # which orientation do we want?
  orient = orient0-orientation
  orient[orient>pi] = orient[orient>pi]-2*pi
  orient[orient<(-pi)] = orient[orient<(-pi)]+2*pi
  
  # define the grating in orientation-frequency space           
  grating = 2* exp(-(log2(f)-log2(freq))^2/bw[1]^2-((orient^2)/bw[2]^2)) # bw is two SDs
  
  # create the Gabor in complex numbers via inverse FFT
  complexGabor = gsignal::fftshift(fft(gsignal::ifftshift(grating, MARGIN = c(1,2)), 
                                       inverse = TRUE), 
                                   MARGIN = c(1,2))
  # get Gabor:
  Gabor = 1 / max(abs(complexGabor)) * (cos(phase)*Real(complexGabor)+sin(phase)*Imag(complexGabor))
  # ... and its dimnames
  dimnames(Gabor) <- list(as.character(y_spat), as.character(x_spat))
  
  return(Gabor)
  
}


# THE ORIGINAL FUNCTION:
# % function image = image_logGabor(imSize,degSize,freq,contrast,center,orientation,L,bw,phase)
# % this function creates a Gabor patch with given orientation, frequency and
# % size (rad and deg of visual angle) in the center of the image
# %
# % 0 phase is the odd positive phase
# %
# % numerical normalization to 1 at peak of the Gabor
# 
# 
# global useGPU
# if isempty(useGPU)
# useGPU = false;
# end
# 
# if ~exist('antiAliasing','var') || isempty(antiAliasing)
# antiAliasing = 4;
# end
# imSize = antiAliasing*imSize;
# 
# 
# if ~exist('f','var') || isempty(f) || ~exist('orient0','var') || isempty(orient0)
# if useGPU
# x = gpuArray.linspace(-imSize(2)/2/degSize(2),(imSize(2)-2)/2/degSize(2),imSize(2));
# y = gpuArray.linspace(-imSize(1)/2/degSize(1),(imSize(1)-2)/2/degSize(1),imSize(1));
# else
#   x = 1:imSize(2);
# y = 1:imSize(1);
# x = (x-ceil(mean(x)))./degSize(2);
# y = (y-ceil(mean(y)))./degSize(1);
# end
# f = sqrt(bsxfun(@plus,x.^2,y'.^2));
#     orient0 = bsxfun(@(x,y) atan2(y,x),x,y');
#          end
#          
#          if useGPU
#          pie = gpuArray(pi);
#          orient = orient0-orientation+pie;
#          %orient(orient>pie)  = orient(orient>pie)-2*pie;
#          %orient(orient<-pie) = orient(orient<-pie)+2*pie;
#          orient = mod(orient,2*pie)-pie;
#          
#          grating = 2* exp(-(log2(f)-log2(freq)).^2./bw(1).^2-(orient.^2./bw(2)^2)); % bw is two standard deviations
#          else
#            orient = orient0-orientation;
#          orient(orient>pi) = orient(orient>pi)-2*pi;
#          orient(orient<-pi) = orient(orient<-pi)+2*pi;
#          
#          grating = 2* exp(-(log2(f)-log2(freq)).^2./bw(1).^2-(orient.^2./bw(2)^2)); % bw is two standard deviations
#          end
#          
#          complexGabor = fftshift(ifft2(ifftshift(grating)));
#          Gabor = 1./max(abs(complexGabor(:))).*(cos(phase).*real(complexGabor)+sin(phase).*imag(complexGabor));
#          
#          center0 = degSize/2;
#          shiftDeg = center-center0;
#          shiftPix = imSize./degSize.*shiftDeg;
#          if shiftPix(1) > 0
#          Gabor((round(shiftPix(1))+1):end,:) = Gabor(1:(end-(round(shiftPix(1)))),:);
#          elseif shiftPix(1)<0
#          Gabor(1:(end+(round(shiftPix(1)))),:) = Gabor((round(-shiftPix(1))+1):end,:);
#          end
#          if shiftPix(2) > 0
#          Gabor(:,(round(shiftPix(2))+1):end) = Gabor(:,1:(end-(round(shiftPix(2)))));
#          elseif shiftPix(2)<0
#          Gabor(:,1:(end+(round(shiftPix(2))))) = Gabor(:,(round(-shiftPix(2))+1):end);
#          end
#          
#          image = L + contrast.*Gabor.*L;
#          
#          image = imresize(image,1/antiAliasing);
#          