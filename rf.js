"use strict"; 

/** 
 *
 * rf.js 
 * Signal processing routines using an emscriptened KissFFT and ROOTjs TGraph's
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 * @license magnet:?xt=urn:btih:1f739d935676111cfff4b4693e3816e664797050&dn=gpl-3.0.txt GPL-v3-or-Later
 * 
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *
 * Sorry, I don't really know javascript, so this is probably all terrible. 
 * Because this relies on KissFFT and jsroot, you'll need to include FFT.js and KissFFT.js ahead of this, e.g:
 *
 *
    <script type="text/javascript" src="jsroot/scripts/JSRootCore.js?2d&hierarchy&io&tree"> </script>
    <script type="text/javascript" src="rfjs/KissFFT.js"></script>
    <script type="text/javascript" src="rfjs/FFT.js"></script>
    <script type="text/javascript" src="rfjs/rf.js"></script>
    <script type="text/javascript" src="myawesomestuff.js"></script>

 */ 



var RF = {}; 
RF.ffts = {}; 

/** Obtain an FFT of the right size 
 *
 * Then you an call forward or inverse on the returned object. 
 *
 **/ 
RF.getFFT =function (size) 
{
  if (this.ffts[size] === undefined)
  {
    this.ffts[size] = new FFTR(size); 
  }

  return this.ffts[size]; 
};


/** performs a forward FFT
 * 
 **/ 
RF.doFFT = function (y) 
{
    var fft = RF.getFFT(y.length); 
    return fft.forward(y).slice(0); 
};


/** Performs inverse FFT 
 *
 * One can optionally specify the number of points in the time domain otherwise will be even. 
 **/ 
RF.doInvFFT = function(Y, Nt = 0) 
{
  var fft = RF.getFFT(Nt == 0 ? Y.length-2 : Nt); 
  return fft.inverse(Y).slice(0); 
};


/** Returns the maximum time and value  of a TGraph
 *
 * Unsigned will use absolute value. 
 *
 * If max_dt is non-zero, will only consider up to max_dt
 *
 * */ 
RF.getMaximumTimeAndValue = function(g, unsigned=true, max_dt = 0) 
{
  var max = 0; 
  var max_t = -1;

  for (var i = 0; i < g.fNpoints; i++) 
  {
    if (max_dt > 0 && Math.abs(g.fX[i]) > max_dt) continue; 
    var val = g.fY[i]; 
    if (unsigned && val < 0) val = -val; 
    
    if (val > Math.abs(max)) 
    {
      max = g.fY[i]; 
      max_t = g.fX[i]; 
    }
  }

  return [max_t, max]; 
}

/** Generates a power spectrum of the given TGraph, returning another TGraph.
 * If you already have the FFT available, you can pass it as Y
 * */ 

RF.makePowerSpectrum = function(g, Y = null) 
{


  var N = g.fX.length;
  if (Y == null) 
  {
    var fft = RF.getFFT(N); 
    Y = fft.forward(g.fY).slice(0); 
  }
  var dt = g.fX[1] - g.fX[0]; 
  var df = 1./(N * dt); 
  var f = []; 
  var P = []; 

  for (var i = 0; i < N/2+1; i++)
  {
    f.push(i*df); 
    var p = (Y[2*i]*Y[2*i] + Y[2*i+1]*Y[2*i+1]) / N; 
    if (i > 0 || i < N/2) p*=2; 
    P.push(10*Math.log10(p)); 
  }

  var G = JSROOT.CreateTGraph(N/2+1, f,P); 
  G.fName = g.fName + "_power"; 
  G.fTitle = g.fTitle.substring(g.fTitle.indexOf(',')+1); 
 
  G.InvertBit(JSROOT.BIT(18)); //make it uneditable

  return G;
}




/** upsamples (in place) a TGraph in the fourier domain 
 *
 *  g: The TGraph
 *  factor: the upsample factor (>1) 
 *  Y: If you already have the fourier trans
 *
 *  returns the fourier transform of the new graph (in case you want to do more with it) 
 *
 * */ 
RF.upsample = function (g, factor=2, Y = null) 
{
  if (factor <=1) return RF.doFFT(g.fY); 

  var N = g.fNpoints; 
  var t0 = g.fX[0]; 
  var dt = g.fX[1] - g.fX[0]; 
  if (Y == null) 
  {
    Y = RF.doFFT(g.fY); 
  }

  var newY = new Float32Array( 2*(g.fNpoints * factor / 2 + 1)); 

  for (var i = 0; i < Y.length; i++) 
  {
    newY[i] = Y[i]/N; 
  }

  g.fNpoints*= factor; 
  g.fY =RF.doInvFFT(newY); 

  for (var i = 0; i < g.fNpoints; i++ )
  {
    g.fX[i] = dt / factor * i + t0 ;
  }
  
  return newY;
}


/** Returns the hilbert transform of a TGraph */ 
RF.hilbertTransform = function(g, Y = null) 
{

  if (Y == null) Y = RF.doFFT(g.fY); 
  var Yp = new Float32Array(Y.length); 

  var N = Y.length; 
  for (var i = 0; i < N/2; i++) 
  {

    Yp[2*i] = -Y[2*i+1] / N; 
    Yp[2*i+1] = Y[2*i] / N; 
  }

  var yp = RF.doInvFFT(Yp, g.fNpoints); 
  var gh= JSROOT.CreateTGraph(g.fNpoints, g.fX, yp); 
  gh.fName = g.fName + "_hilbert"; 
  gh.fTitle = "Hilbert Transform of " + g.fTitle; 
  return gh; 
}


/** Creates the hilbert envelope of a TGraph */ 

RF.hilbertEnvelope = function(g, Y = null, gh = null, H = null) 
{
  if (Y == null) Y = RF.doFFT(g.fY); 

  if (gh == null) 
  {
    gh = RF.hilbertTransform(g,Y); 
  }

  if (H == null) 
  {
    H = JSROOT.CreateTGraph(g.fNpoints, g.fX, g.fY); 
  }
  else
  {
    H.fNpoints = g.fNpoints; 
    H.fX = g.fX; 
    H.fY = []; 
  }

  H.fName = g.fName+"_envelope"; 
  H.fTitle = "Hilbert envelope of " + g.fTitle;
  for (var i = 0; i < g.fNpoints; i++)
  {
    H.fY[i] = Math.sqrt( g.fY[i] * g.fY[i] + gh.fY[i] * gh.fY[i]); 
  }

  H.InvertBit(JSROOT.BIT(18)); //make it uneditable
  return H ; 
}



RF.getGraphMean = function(g) 
{
  var sum = 0; 
  for (var i = 0; i < g.fNpoints; i++) sum += g.fY[i]; 
  return sum / g.fNpoints; 
}

RF.getMean = function(y) 
{
  var sum = 0.; 
  for (var i = 0; i < y.length; i++) sum += y[i]; 
  return sum / y.length; 
}

RF.getRMS = function(g) 
{
  var sum2 = 0; 
  var sum = 0; 
  var N = g.fNpoints;

  for (var i = 0; i < g.fNpoints; i++) 
  {
    sum += g.fY[i];
    sum2 += g.fY[i] * g.fY[i];
  }
  var mean = sum/N;

  return Math.sqrt(sum2/N-mean*mean);
}




RF.rectifyGraph = function(g) 
{
  var mean = RF.getGraphMean(g); 
  for (var i = 0; i < g.fNpoints; i++) g.fY[i] -= mean; 

}

RF.rectify= function(y) 
{
  var mean = RF.getMean(y); 
  for (var i = 0; i < y.length; i++) y[i]-=mean; 

}





RF.range = function(start, N, step = 1) 
{
  var ans = new Float32Array(N); 
  for (var i = 0; i < N; i++) 
  {
    ans[i] = start+step*i; 
  }

  return ans; 
}


/* Calculates the cross-correlation of two graphs, returning the cross-correlation the time domain 
 * Assumes they have the same sampling
 * */ 

RF.crossCorrelation = function ( g1, g2, pad=true, upsample_factor = 4, scale = null, min_freq = 0, max_freq = 0) 
{
  
  if (scale == null) scale = RF.getRMS(g1) * RF.getRMS(g2) ;


  /** most likely we need to pad by a factor of 2 */ 

  var u1 = RF.getGraphMean(g1);
  var u2 = RF.getGraphMean(g2);
  var N = Math.max(g1.fNpoints, g2.fNpoints); 
  var dt = g1.fX[1] - g1.fX[0]; 
  if (pad) N*=2; 
  
  scale *= N/2 * N; 
  var df = 1./(N*dt); 

  var y1 = new Float32Array(N); 
  var y2 = new Float32Array(N); 

  for (var i = 0; i < g1.fNpoints;i++) y1[i] = g1.fY[i]-u1; 
  for (var i = 0; i < g2.fNpoints;i++) y2[i] = g2.fY[i]-u2; 


  var Y1 = RF.doFFT(y1); 
  var Y2 = RF.doFFT(y2); 

  var Y = new Float32Array(upsample_factor * Y1.length); 

  for (var i = 0; i < N/2+1;i++) 
  {
    if (min_freq>0 && i*df < min_freq) continue; 
    if (max_freq>0 && i*df > max_freq) break; 

    var re1 = Y1[2*i];
    var re2 = Y2[2*i];
    var im1 = Y1[2*i+1];
    var im2 = Y2[2*i+1]; 

    Y[2*i] = (re1*re2 + im1*im2)/scale;
    Y[2*i+1] = (im1*re2-re1*im2)/scale;
  }

  var y = RF.doInvFFT(Y); 
  var yrotated = new Float32Array (y.length); 
  for (var i = 0; i < y.length; i++) 
  {
    yrotated[i] = y[(i+y.length/2) % y.length];
  }

  N = y.length; 
  dt = dt/upsample_factor; 
  var x = RF.range(-N/2*dt+g1.fX[0]-g2.fX[0], N, dt); 
  var g = JSROOT.CreateTGraph(y.length, x, yrotated); 

  g.InvertBit(JSROOT.BIT(18)); //make it uneditable


  return g; 
}


/** Linearly interpolates a graph at a point, assuming it's evenly spaced */ 
RF.evalEven = function (g, t) 
{

  var t0 = g.fX[0]; 
  var dt = g.fX[1] - g.fX[0]; 
  var l = Math.floor((t-t0)/dt);
  var u = l+1; 
  if (l < 0) return 0; 
  if (l >= g.fNpoints) return 0; 

  var f = (t - (t0+dt*l))/dt;  
  return f * g.fY[u] + (1-f) * g.fY[l]; 
}



/** A mapper defines a relationship between antennas and times. 
 *  dims: an array with the names of the dimensions (e.g. ["phi","theta"] ), used also to detemrine the length (max 2)! 
 *  nants: The number of antennas
 *  computeDeltaTs:  function with prototype (i,j,x,y) which determines the delta t in terms 
 *  of the parameters.  Dimensions higher than used will be passed 0 
 *  usePair: a function with prototype (i,j) that determines if the pair of antennas can be used
 *  canUse: a function with prototype (i,x,y) that determines if the antenna can be used
 *
 */ 

RF.Mapper = function (dims, nants, computeDeltaTs, usePair = function(i,j) { return true; }, canUse = function(i, x, y) { return true; } ) 
{
  var m = {}; 
  m.ndims = dims.length; 
  m.dims = dims 
  m.nant = nants; 
  m.deltaTs = computeDeltaTs; 
  m.usePair = usePair; 
  m.canUse = canUse; 

  m.createHist = function(nx,xmin,xmax, ny,ymin,ymax)
  {
    var hist = null;
    /* make the histogram */ 
    if (this.ndims == 1) 
    {
      hist = JSROOT.CreateHistogram("TH1F", nx); 
    }
    else if (this.ndims == 2) 
    {
      hist = JSROOT.CreateHistogram("TH2F", nx,ny); 
    }

    hist.fXaxis.fXmin = xmin; 
    hist.fXaxis.fXmax = xmax; 
    hist.fXaxis.fTitle = this.dims[0]; 

    if (this.ndims > 1)
    {
      hist.fYaxis.fXmin = ymin; 
      hist.fYaxis.fXmax = ymax; 
      hist.fYaxis.fTitle = this.dims[1]; 
      hist.fZaxis.fTitle = "Average Cross-Correlation"; 
    }
    else
    {
      hist.fYaxis.fTitle = "Average Cross-Correlation"; 

    }

    return hist;
  }

  m.coherentSum = function(raw_graphs, X, upsample_factor =3, reverse_sign= false) 
  {

    if (raw_graphs.length != this.nant) return null; 

    //check to make sure not all graphs are null! 
    var first_non_null = -1; 
    for (var i = 0; i < raw_graphs.length; i++) 
    {
      if (raw_graphs[i] != null) 
      {
        first_non_null = i; 
        break; 
      }
    }
    if (first_non_null < 0) return null; 


    var graphs = []; 

    // if upsampling, make an upsampled copy 
    //
    if (upsample_factor > 0) 
    {
      for (var i = 0; i < raw_graphs.length; i++) 
      {
        if (raw_graphs[i] == null) graphs.push(null); 
        var gg = JSROOT.CreateTGraph(raw_graphs[i].fNpoints, raw_graphs[i].fX.slice(0), raw_graphs[i].fY.slice(0)); 
        RF.upsample(gg, upsample_factor+1); 
        graphs.push(gg); 
      }
    }
    else
    {
      graphs = raw_graphs; 
    }

    //use first nonnull as time reference

    var idx = first_non_null; 

    var min = graphs[idx].fX[0];
    var max = graphs[idx].fX[graphs[idx].fNpoints-1];
    var dt = (max-min)/(graphs[idx].fNpoints -1); 

    var times = []; 
    for (var i = 0; i < graphs.length; i++) 
    {
      if (graphs[i]==null)
      {
        times.push(null); 
      }
      var this_time = (i == idx ? 0 : this.deltaTs(idx,i, X[0], this.ndims > 1 ? X[1] : 0)); 
      if (reverse_sign) this_time *=-1; 
      times.push(this_time); 
//      console.log(i,X,this_time); 

      if (graphs[i].fX[0] -this_time< min) min = graphs[i].fX[0]-this_time;
      if (graphs[i].fX[graphs[i].fNpoints-1] - this_time > max) max = graphs[i].fX[graphs[i].fNpoints-1] - this_time; 


    }

    var N = Math.floor((max-min)/dt+1); 

    var x = RF.range(min,N,dt); 
    var y= new Float32Array(N); 

    for (var i = 0; i < N; i++) 
    {
      for (var j = 0; j < graphs.length; j++)
      {
        if (graphs[j] == null) continue; 
        var val =  RF.evalEven(graphs[j], x[i]-times[j]); 
        if (!isNaN(val)) y[i] +=val; 
      }
    }

    var ans = JSROOT.CreateTGraph(N,x,y); 
    ans.fTitle = "Coherent Sum";
    ans.InvertBit(JSROOT.BIT(18)); 
    return ans; 
  }

  return m; 
}




RF.dotProduct = function(x,y, Nmax = 0) 
{
  if (Nmax == 0) Nmax = x.length; 
  var sum = 0; 
  for (var i = 0; i < Nmax; i++) sum+= x[i] * y[i]; 
  return sum; 
}

RF.crossProduct = function(u,v) 
{

  return [ u[1]*v[2] - u[2]*v[1] ,
           u[2]*v[0] - u[0]*v[2] , 
           u[0] *v[1] - u[1]*v[0] ]; 
}


RF.magnitude = function(x, Nmax) 
{
  return Math.sqrt(RF.dotProduct(x,x,Nmax)); 
}

/* Returns the angle between two vectors, in radians */ 
RF.angleBetween = function(x,y) 
{
  return Math.acos( RF.dotProduct(x,y) / (RF.magnitude(x) * RF.magnitude(y))); 
}


/** Generate N random gaussians */
RF.randomGaussian = function(N, mean = 0,sigma = 1)
{
  var out = new Float32Array(N); 

  //use marsaglia method
  for (var i = 0; i < N/2; i++) 
  {
    var U = 2*Math.random()-1;
    var V = 2*Math.random()-1;
    var S = U*U+V*V; 
    if ( S > 1)
    { 
      i--; 
      continue; 
    }

    var term = Math.sqrt(-2*Math.log(S) / S); 
    out[2*i]= mean + sigma *U * term;

    if (2*i+1 < N) 
    {
      out[2*i+1]= mean + sigma * V * term; 
    }



  }

  return out; 

}

/* From  https://stackoverflow.com/questions/966225/how-can-i-create-a-two-dimensional-array-in-javascript */ 
RF.createArray = function(length) {
    var arr = new Array(length || 0),
        i = length;
    if (arguments.length > 1) {
        var args = Array.prototype.slice.call(arguments, 1);
        while(i--) arr[length-1 - i] = RF.createArray.apply(this, args);
    }
    return arr;
}



RF.wrap = function(x, period = 360, center = 0)
{
  return x - period * Math.floor((x-center+period/2)/period); 
}


RF.XYMapper = function( xs, ys, c=1) 
{

  return RF.mapper ( ["x (ft)","y (ft)"], xs.length,
      function(i,j, x,y) 
      {

        var dist_i = Math.sqrt(Math.pow(xs[i]-x,2) + Math.pow(ys[i]-y,2)); 
        var dist_j = Math.sqrt(Math.pow(xs[j]-x,2) + Math.pow(ys[j]-y,2)); 
        return (dist_i-dist_j)/c; 
      }); 
}


/* One-dimensional angle mapper, 
 *
 *  Takes only one dimension (not antennas) since not defined if
 *  antennas not all on one axis! 
 **/ 

RF.ElevationMapper = function(zs, c =0.3)
{

  return RF.Mapper( ["theta(deg)"], zs.length, 
    function(i,j, theta_deg, unused) 
    {
      var theta = Math.PI / 180. * theta_deg; 
      var dt= (zs[i]-zs[j])*Math.sin(theta)/c; 
//      console.log(theta_deg,i,j, dt); 
      return -dt; 
    } 
  ); 
}


/** Used by angle mapper */ 
RF.Antenna = function(x,y,z, dx=1,dy=0,dz=0, max_phi=180,max_theta=90) 
{
  var ant = {};
  ant.pos = [x,y,z];
  var bore_mag = Math.sqrt(dx*dx+dy*dy+dz*dz);
  ant.bore = [dx/bore_mag,dy/bore_mag,dz/bore_mag]; 
  ant.phi_width = max_phi;
  ant.theta_width = max_theta; 
  return ant; 
}


/* Two-dimensional angle mapper */ 
RF.AngleMapper = function ( ants,  c = 0.3, phi_0 =0, theta_0 = 0 )
{
  var n = ants.length; 
  var can_use = RF.createArray(n,n); 
  var min_phi =  RF.createArray(n); 
  var max_phi = RF.createArray(n); 
  var min_theta = RF.createArray(n); 
  var max_theta = RF.createArray(n);

  for (var i= 0; i< n; i++)
  {
    var bore_phi = 180./Math.PI* Math.atan2(ants[i].bore[1], ants[i].bore[0]) + phi_0; //check the sign
    var bore_theta = 180./Math.PI* Math.asin(ants[i].bore[2]) + theta_0; //check the sign

    if (ants[i].phi_width >=180) 
    {
      min_phi[i] =-180; 
      max_phi[i] = 180 
    }
    else
    {
      min_phi[i] = RF.wrap(bore_phi - ants[i].phi_width,360,0); 
      max_phi[i] = RF.wrap(bore_phi + ants[i].phi_width,360,0); 
    }

    if (ants[i].theta_width >= 90)
    {

      min_theta[i] =-90; 
      max_theta[i] = 90;
 
    }
    else
    {
      min_theta[i] = RF.wrap(bore_theta - ants[i].theta_width,180,0); 
      max_theta[i] = RF.wrap(bore_theta + ants[i].theta_width,180,0); 
    }

    can_use[i][i] = false; 

    for (var j = i+1; j < n; j++) 
    {
      var bore_xy_i = [ants[i].bore[0], ants[i].bore[1], 0] ; 
      var bore_xy_j = [ants[j].bore[0], ants[j].bore[1], 0] ; 
      var bore_rz_i = [RF.magnitude(ants[i].bore,2),0, ants[i].bore[2]] ; 
      var bore_rz_j = [RF.magnitude(ants[j].bore,2),0, ants[j].bore[2]] ; 
      var boresight_phi_angle =  RF.angleBetween(bore_xy_i, bore_xy_j) * 180/Math.PI; 
      var boresight_theta_angle =  RF.angleBetween(bore_rz_i, bore_rz_j) * 180/Math.PI; 

      can_use[i][j] = boresight_phi_angle < ants[i].phi_width && boresight_phi_angle < ants[j].phi_width && boresight_theta_angle < ants[i].theta_width && boresight_theta_angle < ants[j].theta_width; 
      can_use[j][i] = can_use[i][j]; 
//      console.log(i,j,can_use[i][j]); 
    }
  }
   

  return RF.Mapper( ["phi (deg)","theta (deg)"],ants.length , 

    function (i,j, phi_deg, theta_deg)  //deltaTs
    {
      var phi = Math.PI / 180. * (phi_deg - phi_0); 
      var theta = Math.PI / 180. * (theta_deg - theta_0); 
      var dir = [ Math.cos(phi) * Math.cos(theta), Math.sin(phi) * Math.cos(theta), Math.sin(theta) ]; 
      var diff = [ ants[i].pos[0] - ants[j].pos[0], ants[i].pos[1] - ants[j].pos[1], ants[i].pos[2] - ants[j].pos[2] ];
      return -RF.dotProduct(diff,dir) / c; 
    }, 

    function( i, j)  //usePair
    {
      return can_use[i][j]; 

    }, 

    function(i,phi_deg,theta_deg)  //canUse
    {

       if (theta_deg > max_theta[i]  || theta_deg < min_theta[i]) return false; 
       var between = phi_deg <= max_phi[i] && phi_deg >= min_phi[i]; 
       return  (between && min_phi[i] < max_phi[i]) || (!between && min_phi[i] > max_phi[i]);
    }
  );

}


/** Make a cropped copy of g */
RF.cropWave = function(g, tmin, tmax) 
{

  if (g == null) return null; 
  // find tmin, and tmax

  var imin = 0; 
  var imax = g.fNPoints; 

  for (var i = 0; i < g.fNpoints; i++) 
  {
    if (g.fX[i] < tmin) 
    {
      imin = i+1; 
    }
    if (g.fX[i] > tmax)
    {
      imax =  i-1; 
      break; 
    }
  }

  if (imin >= g.fNPoints) return null; 
  if (imax < 0) return null;

//  console.log(tmin,tmax,imin,imax);
  var cropped= JSROOT.CreateTGraph(imax-imin+1, g.fX.slice(imin,imax+1), g.fY.slice(imin,imax+1)); 

  return cropped; 
}




RF.Spectrogram = function(title, ntime, tmin, tmax, nfreq, fmin, fmax) 
{
  this.hist = JSROOT.CreateHistogram("TH2D", ntime, nfreq); 
  this.sum = JSROOT.CreateHistogram("TH2D", ntime, nfreq); 
  this.norm = JSROOT.CreateHistogram("TH1F", ntime); 
  this.tmin = tmin; 
  this.dt = (this.tmax-this.tmin)/(this.ntime-1); 
  this.nt = ntime; 
  this.nf = nfreq;
  this.df = (fmax-fmin)/(nfreq);
  this.fmin  = fmin; 

  this.norm.fXaxis.fXmin = tmin; 
  this.norm.fXaxis.fXmax = tmax; 
  this.norm.fXaxis.fTimeDisplay = 1; 

  this.hist.fTitle = title; 
  this.hist.fXaxis.fXmin = tmin; 
  this.hist.fXaxis.fXmax = tmax; 
  this.hist.fXaxis.fTimeDisplay = 1; 
  this.hist.fYaxis.fXmin = fmin; 
  this.hist.fYaxis.fXmax = fmax; 

  this.sum.fTitle = title; 
  this.sum.fXaxis.fXmin = tmin; 
  this.sum.fXaxis.fXmax = tmax; 
  this.sum.fXaxis.fTimeDisplay = 1; 
  this.sum.fYaxis.fXmin = fmin; 
  this.sum.fYaxis.fXmax = fmax; 






  this.addY = function(y,t) 
  {
    var Y = RF.doFFT(y); 
    var N = y.length; 
    var P = new Float32Array(N/2+1); 

    for (var i = 0; i <N/2+1; i++)
    {
      P[i] = (Y[2*i]*Y[2*i] + Y[2*i+1]*Y[2*i+1]) / N; 
      if (i > 0 || i < N/2) P[i] *=2; 
    }


    for (var i = 0; i < this.nf; i++) 
    {
      var f = this.fmin + (i+0.5) * this.df; 
      this.sum.Fill(t,f, P[i]); 
    }
    this.norm.Fill(t); 
  }


  this.finalize = function() 
  {

    for (var i = 0; i < this.nt; i++) 
    {
      for (var j = 0; j < this.nf; j++) 
      {
        var ibin = (this.nt+2) * (j+1) + i+1;
//        console.log(i,j,ibin, this.hist.getBinContent(ibin), this.norm.getBinContent(i+1)); 
//        this.hist.setBinContent(ibin, 10 * Math.log10(this.hist.getBinContent(ibin)/ this.norm.getBinContent(i+1))); 
        var n = this.norm.fArray[i+1]; 
        var s = this.sum.fArray[ibin]; 
        var r = s/n; 
        var dbish = n  == 0 ? -20 : 10*Math.log10(r); 
        if (dbish < -20) dbish = -20; 
        this.hist.setBinContent(ibin,  dbish); 
      }
    }
    this.hist.fXaxis.fTitle = "readout time"; 
    this.hist.fYaxis.fTitle = "freq (GHz)"; 
    this.hist.fZaxis.fTitle = "dBish"; 
  }

}



/** This sets up an interferometric map with times defined by the mapper. The delta t's and antennas used
 *  are calculated once at the beginning, which is efficient when one keeps
 *  reusing the same map. 
 *
 *   
 *
 * */ 
RF.InterferometricMap = function ( mapper, nx, xmin, xmax, ny=0, ymin=0,ymax=0, xwrap = false, ywrap = false) 
{

  this.mapper = mapper; 
  this.ndims = mapper.ndims; 
 
  this.hist = mapper.createHist(nx,xmin,xmax,ny,ymin,ymax)

  this.soln = RF.createArray(nx,ny); 

  this.usepair = RF.createArray(mapper.nant,mapper.nant);
  this.xwrap = xwrap;
  this.ywrap = ywrap;
  this.channelNames = function(i) { return i.toString(); }; 
  
  this.nx = nx; 
  this.ny = this.ndims  ==1 ? 1 : ny; 
  this.ymin = ymin; 
  this.ymax = ymax; 
  this.xmin = xmin; 
  this.xmax = xmax; 
  this.dx = (xmax-xmin)/nx; 
  this.dy =  ny > 0 ? (ymax-ymin)/ny : 0; 
  this.nant = mapper.nant;
  this.xcorrs = RF.createArray(this.nant, this.nant);
  this.navg = 0; 

  this.restrict_time_range = false; 
  this.tmin = 0;
  this.tmax = 0;

  this.fmin = 0;
  this.fmax = 0;


  this.upsample_factor = 4; 
  this.is_init = false; 


  this.setTimeRange = function(tmin, tmax) 
  {
    this.restrict_time_range = true; 
    this.tmin = tmin; 
    this.tmax = tmax; 
  }

  this.unsetTimeRange = function() 
  {
    this.restrict_time_range = false; 
  }

  this.setFreqRange = function(fmin, fmax) 
  {
    this.fmin = fmin; 
    this.fmax = fmax; 
  }

  this.unsetFreqRange = function() 
  {
    this.fmin=0;
    this.fmax=0;

  }


  this.init = function() 
  {
    if (this.is_init) return; 

    for (var iant = 0; iant < mapper.nant; iant++)
    {
      for (var jant = 0; jant < mapper.nant; jant++)
      {
        this.usepair[iant][jant] = mapper.usePair(iant,jant);
      }
    }


    for (var ix = 0; ix <this.nx; ix++) 
    {
      for (var iy = 0; iy < this.ny; iy++) 
      {
        var x = (ix+0.5)*this.dx + xmin; 
        var y = this.ndims == 1 ? 0 
                : (iy+0.5)*this.dy + ymin; 

        this.soln[ix][iy] = []; 

        for (var iant  = 0; iant < mapper.nant; iant++)
        {
          if (!mapper.canUse(iant,x,y)) continue; 

          for (var jant = iant+1; jant < mapper.nant; jant++) 
          {
            if (!mapper.canUse(jant,x,y)) continue; 

            if (!this.usepair[iant][jant]) continue;  

            
            var this_soln = {}; 
            this_soln.i = iant; 
            this_soln.j = jant; 
            this_soln.dt = mapper.deltaTs(iant,jant,x,y);
            this.soln[ix][iy].push(this_soln); 
          }
        }
      }
    }
    this.is_init = true; 
  }

  this.setTitle = function (title, xtitle=null, ytitle=null) 
  {
    this.hist.fTitle = title;
    if (xtitle !=null) this.hist.fXaxis.fTitle = xtitle;
    if (xtitle !=null) this.hist.fYaxis.fTitle = ytitle;
  }


  this.deltaTHist = function (i,j)
  {
    this.init(); 
   
    if (i > j) 
    {
      var tmp  = i; 
      i = j; 
      j = tmp; 
    }

    var h = this.mapper.createHist(this.nx,this.xmin,this.xmax,this.ny,this.ymin,this.ymax);
    h.fTitle = "#Delta T (" + i + "," + j+")"; 

    for (var ix = 0; ix < this.nx; ix++)
    {
      for (var iy = 0; iy < this.ny; iy++)
      {
        var val = 0; 
        for (var ipair = 0; ipair < this.soln[ix][iy].length; ipair++)
        {
          if (this.soln[ix][iy][ipair].i == i && this.soln[ix][iy][ipair].j == j)
          {
            val = this.soln[ix][iy][ipair].dt; 
          }
        }
         var ibin = ndim == 1 ? ix+1 : ((this.nx+2) * (iy+1) + ix+1);
         h.setBinContent(ibin, val); 
      }
    }

    return h; 
  }

  this.drawXCorrs = function(where, style_fn = null, max_delay=0, zoom=0, h2_style='margin-top:40%;width:100%;text-align:center;') 
  {
    this.init(); 
    var disp = new JSROOT.GridDisplay(where, "grid" + this.nant + "x" + this.nant); 

    for (var i = 0; i < this.nant; i ++) 
    {
      for (var j = i; j < this.nant; j++) 
      {
        var xcorr_frame = where+"_"+(i + this.nant * j).toString(); 

        var g = i == j ?  RF.crossCorrelation(this.channels[i], this.channels[i]) : this.xcorrs[i][j]; 
        if (style_fn!=null) style_fn(g); 
        g.fTitle = this.channelNames(i) + " WITH  " + this.channelNames(j); 
        var histo = JSROOT.CreateHistogram("TH1I",100); 
        histo.fName = g.fName + "_h";
        histo.fTitle = g.fTitle;
        histo.fXaxis.fXmin = g.fX[0]; 
        histo.fXaxis.fXmax = g.fX[g.fNpoints-1]; 
        histo.fYaxis.fXmin = -1.1;
        histo.fYaxis.fXmax = 1.1;
        histo.fMinimum = -1.1;
        histo.fMaximum = 1.1;
        histo.fXaxis.fTitle = "ns"; 
        histo.fYaxis.fTitle = "correlation"; 
          
        g.fHistogram = histo; 



        JSROOT.draw(xcorr_frame,g, "al", function (painter) 
            {
               var hist = painter.GetObject().fHistogram; 
               painter.root_pad().fGridx = 1; 
               painter.root_pad().fGridy = 1; 
               if (zoom)
               {
                 painter.frame_painter().Zoom("x",-zoom,zoom); 
               }
               else if (max_delay) 
               {
                 painter.frame_painter().Zoom("x",-max_delay,max_delay); 
               }
               JSROOT.redraw(painter.divid, hist, ""); 
            }
            
            ); 

        if (i!=j) 
        {
          //find the maximum or minimum 
          var vals = RF.getMaximumTimeAndValue(this.xcorrs[i][j], true, max_delay); 

          var color_frame = document.getElementById(where+"_"+(j + this.nant*i).toString()); 
          color_frame.innerHTML = "<p> "+this.channelNames(j)+" WITH "+ this.channelNames(i) + " </p><h2 style='"+h2_style+"''> corr<sub>max</sub>="+vals[1].toFixed(4)+"<br>t= " + vals[0].toFixed(2)+ "</h2>"; 
          var deg = (255-Math.floor(Math.abs(vals[1])*255)).toString(16); 
          if (deg.length<2 ) deg = "0"+deg; 
          var string = vals[1] > 0 ? "#" + "ff" + deg+deg : " #" + deg +deg+"ff"; 
//          console.log(string); 
          color_frame.style.backgroundColor = string; 


        }

      }
    }

  }

  this.compute = function(channels, avg = false, reverse_sign = false) 
  {
    this.init(); 

    this.xcorrs = RF.createArray(this.nant, this.nant); 

    var cropped_channels = []; 

    if (this.restrict_time_range) 
    {
      for (var i = 0; i < channels.length; i++) 
      {
        cropped_channels.push(RF.cropWave(channels[i], this.tmin, this.tmax)); 
      }
      channels = cropped_channels;
    }

    this.channels = channels; 

    for (var iant = 0; iant < this.nant; iant++) 
    {
      if (channels[iant] == null) continue;
      for (var jant = iant+1; jant < this.nant; jant++) 
      {
        if (channels[jant] == null) continue;
        if (this.usepair[iant][jant]) 
        {
          var G = RF.crossCorrelation(channels[iant], channels[jant],true,this.upsample_factor,null,this.fmin,this.fmax); 
          this.xcorrs[iant][jant] = G; 
        }
      }
    }

    for (var ix = 0; ix < this.nx; ix++) 
    {
      for (var iy = 0; iy < this.ny; iy++) 
      {
        var sum = 0; 
        var N = this.soln[ix][iy].length; 

        var norm = N; 
        for (var ipair = 0; ipair < N; ipair++) 
        {
          var soln = this.soln[ix][iy][ipair]; 
          if (channels[soln.i] == null || channels[soln.j] == null) 
          {
            norm--; 
            continue; 
          }
          var this_dt = soln.dt;
          sum += RF.evalEven(this.xcorrs[soln.i][soln.j], reverse_sign ? -this_dt: this_dt); 
        }

//        console.log(ix,iy, sum,norm);
        var ibin = (this.ndims == 1 ? ix +1 : (this.nx+2) * (iy+1) + ix+1) ;

        var add_to = avg ? (this.ndims == 1 ? this.hist.getBinContent(ix+1) : this.hist.getBinContent(ix+1,iy+1)) * this.navg : 0;

        var nsum = avg ? (this.navg+1) : 1; 
        var val = norm ? sum/norm : 0; 
//        console.log(ix, iy, ibin,val); 
        this.hist.setBinContent(ibin, (add_to + val)/nsum) ;
      }
    }

    if (avg) this.navg++; 
    else this.navg = 0; 
  }

  this.getMaxes = function(N = 1, min_distance = 10) 
  {
    var maxes = []; 
    var hist = this.hist; 
    for (var imax = 0; imax < N; imax++) 
    {
      var max = 0; 
      var max_x = 0; 
      var max_y = 0; 

      for (var i = 1; i<= hist.fXaxis.fNbins; i++) 
      {
        for (var j = 1; j <= (this.ndims > 1 ? hist.fYaxis.fNbins : 1); j++) 
        {

           var val = this.ndims > 1 ? hist.getBinContent(i,j) : hist.getBinContent(i); 

           if (val > max) 
           {
              var x = hist.fXaxis.GetBinCenter(i);
              var y = this.ndims > 1 ?  hist.fYaxis.GetBinCenter(j) : 0;
              if (maxes.length > 0) 
              {
                var too_close = false; 

                //check that we're not too close to a max

                for (var iimax = 0; iimax< maxes.length; iimax++) 
                {
                  var xdiff = maxes[iimax].x-x; 

                  if (this.xwrap) 
                  {
                    if (xdiff > (this.xmax-this.xmin)/2) 
                    {
                      xdiff = (this.xmax-this.xmin)-xdiff;
                    }

                    if (xdiff < -(this.xmax-this.xmin)/2) 
                    {
                      xdiff = (this.xmax-this.xmin)+xdiff;
                    }

                  }

                  var ydiff = this.ndims > 1 ?  maxes[iimax].y-y : 0; 

                  if (this.ndims > 1 && this.ywrap) 
                  {
                    if (ydiff > (this.ymax-this.ymin)/2) 
                    {
                      ydiff = (this.ymax-this.ymin)-ydiff;
                    }

                    if (ydiff < -(this.ymax-this.ymin)/2) 
                    {
                      ydiff = (this.ymax-this.ymin)+ydiff;
                    }

                  }

                  if (xdiff*xdiff + ydiff*ydiff < min_distance*min_distance)
                  {
                    too_close = true; 
                    break; 
                  }
                }

                if (too_close) continue; 
              }

              max = val; 
              max_x = x; 
              max_y = y; 
            }
        }

      }

//      console.log(max_x,max_y,max); 

      var this_max= { x: max_x, y: max_y, max: max};
      maxes.push(this_max); 

    }
    return maxes; 
  }

}


/** Perform an IIR filter (in-place) on a TGraph (g) 
 * using coefficients b and a
 * */
RF.IIRFilter = function(g, b,a) 
{

  if (a == null || a.length == 0) a = [1]; 
  if (b == null || b.length == 0) b = [1]; 

  var yNew = new Float32Array(g.fNpoints); 

  var inv = 1./a[0]; 
  for (var i = 0; i < g.fNpoints; i++) 
  {
    var val = 0; 

    for (var j = 0; j < Math.min(i+1,b.length); j++)
    {
      val += b[j] * g.fY[i-j]; 
    }

    for (var k = 1; k < Math.min(i+1,a.length); k++) 
    {
      val -= a[k] * yNew[i-k]; 
    }

    yNew[i] = val *inv; 

  }

  for (var i = 0; i < g.fNpoints;i++) 
  {
    g.fY[i] = yNew[i]; 
  }
}




RF.shiftTimes = function(g, t) 
{
  if (g==null) return; 
  for (var i = 0; i < g.fNpoints; i++) { g.fX[i] +=t; }
}

RF.sumWithDelays = function( raw_graphs, times,upsample = 3 ) 
{

  if (!(raw_graphs.length >0) || raw_graphs.length!=times.length) return null; 


  var graphs = []; 

  // if upsampling, make an upsampled copy 
  //
  if (upsample > 0) 
  {
    for (var i = 0; i < raw_graphs.length; i++) 
    {
      if (raw_graphs[i] == null) graphs.push(null); 
      var gg = JSROOT.CreateTGraph(raw_graphs[i].fNpoints, raw_graphs[i].fX.slice(0), raw_graphs[i].fY.slice(0)); 
      RF.upsample(gg, upsample+1); 
      graphs.push(gg); 
    }
  }
  else
  {
    for (var i = 0; i < raw_graphs.length; i++) graphs.push(raw_graphs[i]); 
  }


  //set up the grid, use min / max and dt of first graph 

  //find the first non-null graph

  var idx = -1;
  for (var i = 0; i < graphs.length; i++) 
  {
    if (graphs[i]!=null)
    {
      idx = i; 
      break; 
    }
  }
  if (idx == -1) 
  {
    return null;// all null
  }

  var min = graphs[idx].fX[0] - times[idx]; 
  var max = graphs[idx].fX[graphs[idx].fNpoints-1] - times[idx]; 
  var dt = (max-min)/(graphs[idx].fNpoints -1); 


  for (var i = 0; i < graphs.length; i++) 
  {
    if (graphs[i]==null) continue;
    if (graphs[i].fX[0] -times[i]< min) min = graphs[i].fX[0]-times[i]; 
    if (graphs[i].fX[graphs[i].fNpoints-1] > max) max = graphs[i].fX[graphs[i].fNpoints-1] - times[i]; 
  }
  var N = Math.floor((max-min)/dt+1); 

  var x = RF.range(min,N,dt); 
  var y= new Float32Array(N); 

  for (var i = 0; i < N; i++) 
  {
    for (var j = 0; j < graphs.length; j++)
    {
      if (graphs[j] == null) continue; 
      var val =  RF.evalEven(graphs[j], x[i]-times[j]); 
      if (!isNaN(val)) y[i] +=val; 
    }
  }

  var ans = JSROOT.CreateTGraph(N,x,y); 
  ans.fTitle = "Coherent Sum";
  ans.InvertBit(JSROOT.BIT(18)); 
  return ans; 

}
