"use strict"; 

/** 
 *  Various things for working with RF TGraphs
 *
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


/** performs FFT */ 

RF.doFFT = function (y) 
{
    var fft = RF.getFFT(y.length); 
    return fft.forward(y); 
};


/** Performs inverse FFT */ 
RF.doInvFFT = function(Y, Nt = 0) 
{
  var fft = RF.getFFT(Nt == 0 ? Y.length : Nt); 
  return fft.inverse(Y); 
};


/** Generates a power spectrum of the given TGraph, returning another TGraph.
 * If you already have the FFT available, you can pass it as Y
 * */ 

RF.makePowerSpectrum = function(g, Y = null) 
{


  var N = g.fX.length;
  if (Y == null) 
  {
    var fft = RF.getFFT(N); 
    Y = fft.forward(g.fY); 
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


/** Returns the hilbert transform of g */ 
RF.hilbertTransform = function(g, Y = null) 
{

  if (Y == null) Y = RF.doFFT(g.fY); 
  var Yp = new Float32Array(Y.length); 

  var N = g.fNpoints; 
  for (var i = 1; i < g.fNpoints/2; i++) 
  {

    Yp[2*i] = -Y[2*i+1] / N; 
    Yp[2*i+1] = Y[2*i] / N; 
  }

  var yp = RF.doInvFFT(Yp); 
  var gh= JSROOT.CreateTGraph(g.fNpoints, g.fX, yp); 
  gh.fName = g.fName + "_hilbert"; 
  gh.fTitle = "Hilbert Transform of " + g.fTitle; 
  return gh; 
}


/** Creates the hilbert envelope of a graph */ 

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


RF.getMean = function(g) 
{
  var sum = 0; 
  for (var i = 0; i < g.fNpoints; i++) sum += g.fY[i]; 
  return sum / g.fNpoints; 
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

  return Math.sqrt(sum2/N-mean);
}



RF.rectify = function(g) 
{
  var mean = RF.getMean(g); 
  for (var i = 0; i < g.fNpoints; i++) g.fY[i] -= mean; 

}


RF.range = function(start, N, step = 1) 
{
  var ans = Array(N); 
  for (var i = 0; i < N; i++) 
  {
    ans[i] = start+step*i; 
  }

  return ans; 
}


/* Calculates the cross-correlation of two graphs, returning the cross-correlation the time domain 
 * Asumes they have the same sampling... 
 * */ 

RF.crossCorrelation = function ( g1, g2, pad=true, upsample = 4) 
{
  
  /** most likely we need to pad by a factor of 2 */ 

  var N = Math.max(g1.fNpoints, g2.fNpoints); 
  var dt = g1.fX[1] - g1.fX[0]; 
  if (pad) N*=2; 
  

  var y1 = new Float32Array(N); 
  var y2 = new Float32Array(N); 

  for (var i = 0; i < g1.fNpoints;i++) y1[i] = g1.fY[i]; 
  for (var i = 0; i < g2.fNpoints;i++) y2[i] = g2.fY[i]; 

  var Y1 = RF.doFFT(y1); 
  var Y2 = RF.doFFT(y2); 
  Y = new Float32Array(upsample * Y1.length); 

  for (var i = 0; i < N/2+1;i++) 
  {
    var re1 = Y1[2*i]; 
    var re2 = Y2[2*i]; 
    var im1 = Y1[2*i+1]; 
    var im2 = Y2[2*i+1]; 

    Y[2*i] = re1*re1 + im1*im2;
    Y[2*i+1] = im1*re2-re1*im2; 
  }

  var y = RF.doInvFFT(Y); 
  y = y.slice(N,y.length).concat(y.slice(0,N)); 
  var x = RF.range(-N*dt, 2*N, dt); 
  var g = JSROOT.CreateTGraph(y.length, x, y); 

  return g; 
}


/** Interpolates a graph at a point, assuming it's evenly spaced */ 
RF.evalEven = function (g, t) 
{

  var t0 = g.fX[0]; 
  var dt = g.fX[1] - g.fX[0]; 
  var l = Math.floor((t-t0)/dt);
  var u = l+1; 
  if (l < 0) return 0; 
  if (l >= g.fNpoints) return 0; 

  var f = t - (t0+dt*l) 
  return f * g.fY[u] + (1-f) * g.fY[l]; 
}



RF.Mapper = function (nants, computeDeltaTs, usePair = function(i,j) { return true; }, canUse = function(i, x, y) { return true; } ) 
{
  var m = {}; 
  m.nants = nants; 
  m.deltaTs = computeDeltaTs; 
  m.usePair = usePair; 
  m.canUse = canUse; 
  return m; 
}




RF.Antenna = function(x,y,z, dx,dy,dz, max_phi,max_theta) 
{
  var ant = {};
  ant.pos = [x,y,z];
  ant.bore = [dx,dy,dz]; 
  ant.phi_width = phi_width;
  ant.theta_width = theta_width; 
  return ant; 
}


RF.dotProduct = function(x,y, Nmax = 0) 
{
  if (Nmax == 0) Nmax = x.length; 
  var sum = 0; 
  for (var i = 0; i < Nmax; i++) sum+= x[i] * y[i]; 
  return sum; 
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

/* From  https://stackoverflow.com/questions/966225/how-can-i-create-a-two-dimensional-array-in-javascript */ 
RF.createArray = function(length) {
    var arr = new Array(length || 0),
        i = length;
    if (arguments.length > 1) {
        var args = Array.prototype.slice.call(arguments, 1);
        while(i--) arr[length-1 - i] = createArray.apply(this, args);
    }
    return arr;
}



RF.wrap = function(x, period = 360, center = 0)
{
  return x - period * Math.floor((x-center+period/2)/period); 
}

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

    min_phi[i] = RF.wrap(bore_phi - ant[i].phi_width,360,180); 
    max_phi[i] = RF.wrap(bore_phi + ant[i].phi_width,360,180); 
    min_theta[i] = RF.wrap(bore_theta - ant[i].phi_width,180,0); 
    max_theta[i] = RF.wrap(bore_theta + ant[i].phi_width,180,0); 

    can_use[i][i] = false; 

    for (var j = i+1; j < n; j++) 
    {
      var bore_xy_i = [ants[i].bore[0], ants[i].bore[1], 0] ; 
      var bore_xy_j = [ants[j].bore[0], ants[j].bore[1], 0] ; 
      var bore_rz_i = [RF.magnitude(ants[i].bore,2),0, ants[i].bore[2]] ; 
      var bore_rz_j = [RF.magnitude(ants[j].bore,2),0, ants[j].bore[2]] ; 
      var boresight_phi_angle =  angleBetween(bore_xy_i, bore_xy_j) * 180/Math.PI; 
      var boresight_theta_angle =  angleBetween(bore_rz_i, bore_rz_j) * 180/Math.PI; 

      can_use[i][j] = boresight_phi_angle < ants[i].phi_width && boresight_phi_angle < ants[j].phi_width && boresight_theta_angle < ants[i].theta_width && boresight_theta_angle < ants[j].theta_width; 
      can_use[j][i] = can_use[i][j]; 
    }
  }
   

  return RF.Mapper( ants.length , 

    function (i,j, phi_deg, theta_deg)  //deltaTs
    {
      var phi = Math.PI / 180. * (phi_deg - phi_0); 
      var theta = Math.PI / 180. * (theta_deg - theta_0); 
      var dir = [ cos(phi) * cos(theta), sin(phi) * cos(theta), sin(theta) ]; 
      var diff = [ ants[i].pos[0] - ants[j].pos[0], ants[i].pos[1] - ants[j].pos[1], ants[i].pos[2] - ants[j].pos[2] ];
      var distance = RF.magnitude(diff); 
      var angle = RF.angleBetween(dir,diff); 
      return c * cos(angle) * distance; 
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




/** This sets up an interferometric map (which will be a TH2) 
 *  with times defined by the mapper. The delta t's and antennas used
 *  are calculated once at the beginning, which is efficient when one keeps
 *  reusing the same map
 *
 * */ 
RF.InterferometricMap = function ( nx, xmin, xmax, ny, ymin,ymax, mapper) 
{

  var map = {}; 

  /* make the histogram */ 
  map.hist = JSROOT.CreateHistogram("TH2F", nx,ny); 
  map.hist.fXaxis.fXmin = xmin; 
  map.hist.fXaxis.fXmax = xmax; 
  map.hist.fXaxis.fYmin = ymin; 
  map.hist.fXaxis.fYmax = ymax; 

  map.soln = RF.createArray(nx,ny); 

  map.usepair = RF.createArray(nx,ny);
  
  map.nx = nx; 
  map.ny = ny; 
  map.ymin = ymin; 
  map.ymax = ymax; 
  map.xmin = xmin; 
  map.xmax = xmax; 
  map.dx = (xmax-xmin)/nx; 
  map.dy = (ymax-ymin)/ny; 
  map.nant = mapper.nant

  for (var iant = 0; iant < mapper.nant; iant++)
  {
    for (var jant = 0; jant < mapper.nant; jant++)
    {
      usepair[iant][jant] = mapper.usePair(iant,jant);
    }
  }


  for (var ix = 0; x < nx; ix++) 
  {
    for (var iy = 0; iy < ny; iy++) 
    {
      var x = ix*map.dx + xmin; 
      var y = iy*map.dy + ymin; 

      map.soln[i][j] = []; 

      for (var iant  = 0; iant < mapper.nant; iant++)
      {
        if (!mapper.canUse(iant,x,y)) continue; 

        for (var jant = iant+1; jant < mapper.nant; jant++) 
        {
          if (!mapper.canUse(jant,x,y)) continue; 
          if (usepair[iant][jant]) continue;  

          var this_soln = {}; 
          this_soln.i = iant; 
          this_soln.j = jant; 
          this_soln.dt = mapper.deltaTs(iant,jant,x,y);
          map.soln[i][j].push(this_soln); 
        }
      }
    }
  }

  map.compute = function(channels) 
  {

    var xcorrs = RF.createArray(this.nants, this.nant); 

    for (var iant = 0; iant < this.nant; iant++) 
    {
      for (var jant = iant+1; jant < this.nant; jant++) 
      {
        if (this.usepair[iant][jant]) 
        {
          xcorrs[iant][jant] = RF.crossCorrelation(channels[iant], channels[jant]); 
        }
      }
    }

    for (var ix = 0; ix < this.nx; ix++) 
    {
      for (var iy = 0; iy < this.ny; iy++) 
      {
        var sum = 0; 
        var norm = this.solns[ix][jx].length; 

        for (var ipair = 0; ipair < norm; ipair++) 
        {
          var soln = this.solns[ix][jx][ipair]; 
          sum += RF.evalEven(xcorrs[soln.i][soln.j], soln.dt); 
        }

        this.hist.setBinContent(ix,iy, norm? sum/norm : 0); 
      }
    }
  }

  return map; 
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

