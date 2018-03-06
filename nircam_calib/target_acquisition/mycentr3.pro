function mycentr3,image,hwm,checkbox_hw,iter,initial_center
;calculates the simple centroid for an image in both the X and Y direction
;Then introduces a weighting scheme for recentering 
;
;INPUTS:
;image - the TA image (2D array)
;hwm  - the half width of the centroid box
;checkbox_hw - half width of the checkbox used for the
;              initial estimate of the peak location.
;              For NIRCam we use a 3x3 box, so set
;              checkbox_hw to 1. (i.e peak-1:peak+1)
;iter - the number of iterations to use
;initial_center  - a 2 element array. It is only
;                  here as an input so that the peak
;                  as determined by the initial checkbox
;                  algorithm can be returned by reference.
;                  Input value can be anything as it
;                  will be ignored.
;
;OUTPUTS:
;cent - 2 element array giving the final centroid results
;initial_center - 2 element array giving the initial centroid
;                 results from the 3x3 checkbox
;
;HISTORY:
;Written by John Stansberry. Modifications for NIRCam
;TA with saturated pixel study by Bryan Hilbert. Jan 2018
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;delta_checkbox is the cadence in pixels
;pixels of how often the checkbox is applied.
;Should be set to 1 for all instruments, so that
;the checkbox is shifted one pixel at a time.
delta_checkbox = 1
iter=iter-1

;older method. From John's original script
;peak=where(image eq max(image[1:N_elements(image(*,0))-2,1:N_elements(image(*,0))-2]))
;peak=array_indices(image,peak)

;Hilbert's interpretation of the first step of GENTALOCATE
;from Goudfrooji's description.
;find peak signal in 3x3 checkbox
peaksig = 0
peakx = 0
peaky = 0
xl = n_elements(image[*,0])
yl = n_elements(image[0,*])
for i=checkbox_hw,xl-checkbox_hw-1 do begin
   for j=checkbox_hw,yl-checkbox_hw-1 do begin
      sig = total(image[i-checkbox_hw:i+checkbox_hw,j-checkbox_hw:j+checkbox_hw])
;for i=1,n_elements(image[*,0])-2,delta_checkbox do begin
;   for j=1,n_elements(image[0,*])-2,delta_checkbox do begin
;      sig = total(image[i-1:i+1,j-1:j+1])
      if sig gt peaksig then begin
         peaksig = sig
         peakx = i
         peaky = j
      endif
   endfor
endfor
peak = [peakx,peaky]
         

cent=fltarr(2)
old_cent=cent
peak=round(peak)
xsum=0.
ysum=0.
sum=0.
weight=0
sumweight=0
print,'Initial peak location:',peak

for i=peak(0)-hwm,peak(0)+hwm do begin 
   for j=peak(1)-hwm,peak(1)+hwm do begin 
      sum=sum+image(i,j)
      xsum=xsum+i*image(i,j)
      ysum=ysum+j*image(i,j)
   end
end
                    
cent(0)=1.*xsum/sum
cent(1)=1.*ysum/sum
orig=cent
old_cent(0)=0
old_cent(1)=0
v=0
for conv=0,iter do begin 

   for i=floor(cent(0)-hwm),ceil(cent(0)+hwm) do begin 
       v=v+1
      for j=floor(cent(1)-hwm),ceil(cent(1)+hwm) do begin 
         xwei=0
         yweight=0
         xoff=i-cent(0)
         yoff=j-cent(1)
         if (abs((xoff)) le hwm) then begin 
            xweight=1
            endif else begin 
               if ((abs((xoff)) ge hwm) and (abs((xoff)) le hwm+1)) then begin 
                  xweight=hwm+1-abs(xoff)
               endif
            endelse
            if (abs(yoff) le hwm) then begin 
               yweight=1
            endif else begin 
              if (abs(yoff) gt hwm) and (abs(yoff) le hwm+1) then begin 
                 yweight=hwm+1-abs(yoff)
              endif
           endelse
           if (abs(xoff) gt hwm+1 or abs(yoff) gt hwm+1) then begin 
               xweight=0 
               yweight=0
           end

            weight=xweight*yweight
            sumweight=sumweight+weight
            sum=image(i,j)*weight+sum
            xsum=i*image(i,j)*weight+xsum
            ysum=j*image(i,j)*weight+ysum
         endfor
     endfor

     old_cent=cent
            cent(0)=xsum/sum
            cent(1)=ysum/sum

        endfor

;function returns the calculated centroid value
;but also passes back the initial peak value by
;reference
initial_center = peak
return,cent

end
