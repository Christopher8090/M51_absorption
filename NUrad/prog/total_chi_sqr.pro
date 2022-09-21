pro total_chi_sqr,model,qyear,scaabs,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,f_uv,f_uv4,f_uv6,f_uv7,f_BVIK,f_BVIK3,f_BVIK5
compile_opt idl2

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd
print_text = 0  ; to print information about each step, set to 1.

root = '/nfs/d58/vlk/sedmodel/cinman/m51a/NUrad_M51a/'
save_dir = root+'saves/model/'

wavelength = STRARR(18)
wavelength[0]  = 'uv15'
wavelength[1]  = 'uv22'
wavelength[2]  = 'uv36'
wavelength[3]  = 'b'
wavelength[4]  = 'v'
wavelength[5]  = 'i'
wavelength[6]  = 'j'
wavelength[7]  = 'k'
wavelength[8]  = '3.6um'
wavelength[9]  = '4.5um'
wavelength[10] = '5.8um'
wavelength[11] = '8.0um'
wavelength[12] = '24um'
wavelength[13] = '70um'
wavelength[14] = '160um'
wavelength[15] = '250um'
wavelength[16] = '350um'
wavelength[17] = '500um'

chi_wave = dblarr(n_elements(wavelength))
chi_n = dblarr(n_elements(wavelength))
n_annuli = 250
read, combined, prompt='combine LMC and MW dust? [1]=yes: '
if combined eq 1 then model = 'combined'

for i = 0, n_elements(wavelength)-1 do begin
	;if i eq 12 then continue
	if i eq 0 then print, 'band    chi^2_r            chi^2          n_annuli'
	fname = save_dir+model+'_datafile_'+wavelength[i]+'_'+scaabs+'.save'
	RESTORE, fname
	chi_wave[i] = total(chi_sqr,/nan)
	chi_n[i] = chi_count
	print, wavelength[i], chi_sqr_r, chi_wave[i], chi_n[i]
endfor

sum_chi_wave = total(chi_wave)
sum_n_annuli = total(chi_n)
;print, 'total:', sum_chi_wave, sum_n_annuli
total_chi_sqr = sum_chi_wave/ sum_n_annuli
print, ''
print, 'total chi^2_r for '+model+' model: ',total_chi_sqr
stop
end
