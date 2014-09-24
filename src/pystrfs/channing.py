import sys
from pystrfs import *

from browser.model.analysis.strflab.linear_rf import LinearRF

from scipy import signal

all_rfs = meta.Session.query(LinearRF).filter_by(preprocess_id=1)

power_func = lambda psth: (psth**2).sum()
snr = 3.0


def get_channing_units():

    ulst = []

    f = open('/auto/fhome/mschachter/channing_units.txt')
    ustrs = f.readlines()
    f.close()
    uids = [int(x.strip()) for x in ustrs if len(x.strip()) > 0]

    for uid in uids:
        q = meta.Session.query(Unit).filter_by(id=uid)
        unit = q.one()
        ulst.append(unit)

    return ulst

def write_channing_units():

    for u in get_channing_units():
        create_channing_unit_file(u)

def get_all_unit_data(unit):

    classresps = [classresp for classresp in unit.class_resps if classresp.class_name == "masked song"][0] #Get masked song class response object

    print "Collecting all responses for unit %d" % unit.id
    #masked_resps = [r for r in classresps.repeated] #Get all masked response objects
    #unmasked_resps = [r.paired for r in masked_resps] #Get corresponding clean song response objects
    #masker_resps = [[s.contained().repeated for s in r.singles if s.contained() is not None] for r in masked_resps] #Get corresponding noise response objects

    unit_dict = dict()

    #get song STRF
    unit_strfs = LinearRF.unit_id.in_([unit.id])
    query = all_rfs.filter(unit_strfs).filter_by(class_name='song')
    strf = query.one()

    # Spectrogram parameters
    print "\tGenerating spectrogram parameters"
    twindow = strf.preprocesser.t_window
    sr = classresps.repeated[0].presentation.stim.file.get_samplerate()
    resampled_sr = np.ceil(np.double(sr) / 1000.0) * 1000.0 #sr is resampled to be a multiple of 1000 when computing the spectrogram
    winlength = np.fix(twindow * resampled_sr)
    winlength = np.fix(winlength / 2.0) * 2.0
    unit_dict["spectrogram_params"] = {"twindow": twindow, "winlength": winlength, "sr": sr}

    for jj, mresp in enumerate(classresps.repeated):

        if mresp.singles[0].contained() is None:
            break
        print "\tResponse #%d" %jj
        # Initializations
        print "\t\tInitializing"
        for s in ["masked", "masker"]:
            for ss in ["resp", "pred_resp", "wav", "spec"]:
                unit_dict[s][ss].append([])
        for ss in ["SNR", "SNR_pred", "start"]:
            unit_dict["masked"][ss].append([])

        uresp = mresp.paired


        # Make PSTHs
        print "\t\tMaking unmasked and masked PSTHs"
        if len(uresp.spiketimes):
            upsth = uresp.make_psth().data
            upower = power_func(upsth)
        else:
            upsth = []
            upower = 0
            print "Unit %d did not spike to unmasked stimulus %d"%(unit.id, uresp.stim_number)
            #exitFlag = True
            #break

        if len(mresp.spiketimes):
            mpsth = mresp.make_psth().data
        else:
            mpsth = []
            print "Unit %d did not spike to masked stimulus %d"%(unit.id, mresp.stim_number)

        unit_dict["masked"]["resp"][jj] = mpsth
        unit_dict["unmasked"]["resp"].append(upsth)



        # Get signal sound data
        print "\t\tGetting unmasked sound data"
        uwav = mresp.presentation.stim.file.get_data()
        unit_dict["unmasked"]["wav"].append(uwav)

        uspec, time, freq = spectrogram(uwav, sr, window_length=winlength)
        pass_freqs = (freq >= 250) & (freq <= 8100)
        uspec = log_spec(uspec[pass_freqs])
        unit_dict["unmasked"]["spec"].append(uspec)

        # Create predicted PSTH
        print "\t\tComputing unmasked predicted PSTH"
        uspec = normalize_spec(uspec) #Need to figure out how normalization is done.
        upsth_pred = conv_activations_2d(uspec, strf.strf.data)[0].squeeze()
        unit_dict["unmasked"]["pred_resp"].append(upsth_pred)

        for kk, sresp in enumerate(mresp.singles):
            print "\t\tResponse single #%d"%kk
            try:
                nresp = sresp.contained().repeated
            except AttributeError:
                break

            # Make PSTH
            print "\t\t\tComputing masker PSTH"
            if len(nresp.spiketimes):
                npsth = nresp.make_psth().data
                npower = power_func(npsth)
            else:
                npsth = []
                npower = 0
                print "Unit %d did not spike to noise stimulus %d"%(unit.id, nresp.stim_number)
            unit_dict["masker"]["resp"][jj].append(npsth)

            # Calc response SNR
            print "\t\t\tCalculating response SNR"
            if npower != 0:
                SNR = upower / npower
            else:
                SNR = np.NaN
            unit_dict["masked"]["SNR"][jj].append(SNR)

            # Get noise and masked sound data
            print "\t\t\tGetting masker sound data and computing masked sound data"
            nwav = nresp.presentation.stim.file.get_data()
            start = float(sresp.start - sresp.contained().start) * sr
            mwav = add_signal_to_noise(uwav, nwav, snr, start)

            unit_dict["masker"]["wav"][jj].append(nwav)
            unit_dict["masked"]["wav"][jj].append(mwav)
            unit_dict["masked"]["start"][jj].append(float(sresp.start - sresp.contained().start))

            nspec = spectrogram(nwav, sr, window_length=winlength)[0]
            nspec = log_spec(nspec[pass_freqs])
            mspec = spectrogram(mwav, sr, window_length=winlength)[0]
            mspec = log_spec(mspec[pass_freqs])

            unit_dict["masked"]["spec"][jj].append(mspec)
            unit_dict["masker"]["spec"][jj].append(nspec)

            # Create predicted PSTHs
            print "\t\t\tCreating predicted PSTHs to masker and masked spectrograms"
            nspec = normalize_spec(nspec)
            mspec = normalize_spec(mspec)

            npsth_pred = conv_activations_2d(nspec, strf.strf.data)[0].squeeze()
            mpsth_pred = conv_activations_2d(mspec, strf.strf.data)[0].squeeze()

            unit_dict["masked"]["pred_resp"][jj].append(mpsth_pred)
            unit_dict["masker"]["pred_resp"][jj].append(npsth_pred)

            #Calc predicted SNR
            print "\t\t\tComputing predicted SNR"
            upower_pred = power_func(upsth_pred)
            npower_pred = power_func(npsth_pred)
            if npower_pred != 0:
                SNR_pred = upower_pred / npower_pred
            else:
                SNR_pred = np.NaN
            unit_dict["masked"]["SNR_pred"][jj].append(SNR_pred)

            sys.stdout.flush()


def create_channing_unit_file(unit, output_dir='/tmp'):

    unitFile = os.path.join(output_dir, 'unit_%d.h5' % unit.id)

    f = h5py.File(unitFile, 'w')
    f.attrs['electrode'] = int(unit.electrode_site_number)
    f.attrs['ldepth'] = 0.0
    f.attrs['rdepth'] = float(unit.recsite.real_depth_um)
    f.attrs['site'] = '%d' % unit.recsite.id
    f.attrs['source_directory'] = 'database'

    print 'Creating unit file for %d...' % unit.id


    con_class_name = 'song'
    mlnoise_class_name = 'ml noise'
    maskedcon_class_name = 'masked song'

    #write conspecific unmaksed songs/responses
    con_class_resps = [x for x in unit.class_resps if x.class_name == con_class_name]
    if len(con_class_resps) == 0:
        print 'No conspecific songs for unit %d' % unit.id
        return
    elif len(con_class_resps) > 1:
        print 'More than one class response for conspecific songs for unit %d' % unit.id
    con_class_resps = con_class_resps[0]

    write_class_resps(f, con_class_resps)

    #write mlnoise
    mlnoise_class_resps = [x for x in unit.class_resps if x.class_name == mlnoise_class_name]
    if len(mlnoise_class_resps) == 0:
        print 'No ml-noise for unit %d' % unit.id
        return
    elif len(mlnoise_class_resps) > 1:
        print 'More than one class response for ml-noise for unit %d' % unit.id
    mlnoise_class_resps = mlnoise_class_resps[0]

    write_class_resps(f, mlnoise_class_resps)

    #write conspecific masked songs/responses
    maskedcon_class_resps = [x for x in unit.class_resps if x.class_name == maskedcon_class_name]
    if len(maskedcon_class_resps) == 0:
        print 'No masked songs for unit %d' % unit.id
        return
    elif len(maskedcon_class_resps) > 1:
        print 'More than one class response for masked songs for unit %d' % unit.id
    maskedcon_class_resps = maskedcon_class_resps[0]

    write_class_resps_masked(f, maskedcon_class_resps)

    print 'Writing unit hdf5 file to %s...' % unitFile
    f.close()


def write_class_resps(f, class_resps):

    stimGroup = class_resps.class_name.replace(' ','_')
    print 'Writing class responses for %s' % stimGroup
    classGroupName = '%s' % stimGroup
    cpGroup = f.create_group(classGroupName)

    for k,repPres in enumerate(class_resps.repeated):

        stim_number = repPres.presentation.stim.number
        stimGroupName = '%d' % stim_number
        print 'Creating stim group: %s' % stimGroupName
        stimGroup = cpGroup.create_group(stimGroupName)

        stim = repPres.presentation.stim

        stimGroup.attrs['original_wavfile'] = str(repPres.stim_filename)
        stimGroup.attrs['stim_class'] = str(repPres.stim_group)
        stimGroup.attrs['stim_duration'] = float(stim.file.duration)
        stimGroup.attrs['stim_md5'] = str(stim.file.MD5)
        stimGroup.attrs['stim_source'] = 'Unknown'
        stimGroup.attrs['stim_source_sex'] = 'Unknown'
        stimGroup.attrs['stim_type'] = str(stim.stim_group)
        stimGroup.attrs['tdt_wavfile'] = str(stim.filename)
        stimGroup.attrs['tstat'] = 0.0
        stimGroup.attrs['pvalue'] = 1.0
        stimGroup.attrs['zscore'] = 0.0

        spikeTrials = repPres.spiketimes

        for trialNum,spikes in enumerate(spikeTrials):
            trialResponseName = '%d' % int(trialNum+1)
            #print '\tCreating trial response with %d spikes: %s' % (len(spikes), trialResponseName)
            trialResponseGroup = stimGroup.create_group(trialResponseName)
            trialResponseGroup.attrs['bg_rate'] = 0.0
            trialResponseGroup.attrs['peri_rate'] = len(spikes) / float(stim.file.duration)
            if len(spikes) > 0:
                trialResponseGroup['spike_times'] = np.array(spikes)
            else:
                trialResponseGroup['spike_times'] = -999

def write_class_resps_masked(f, class_resps):

    stimGroup = class_resps.class_name.replace(' ','_')
    print 'Writing class responses for %s' % stimGroup
    classGroupName = '%s' % stimGroup
    cpGroup = f.create_group(classGroupName)

    for k,repPres in enumerate(class_resps.repeated):

        song_stim_number = repPres.presentation.stim.number

        for m,single in enumerate(repPres.singles):

            masker = single.contained()
            if masker is None:
                #some weird stuff with birds own song
                continue

            noise_stim_number = masker.stim_number

            stimGroupName = '%d_%d' % (song_stim_number, noise_stim_number)
            print 'Creating stim group: %s' % stimGroupName
            stimGroup = cpGroup.create_group(stimGroupName)
            stimGroup.attrs['stim_duration'] = float(masker.duration)
            stimGroup.attrs['song_stim_number'] = song_stim_number
            stimGroup.attrs['noise_stim_number'] = noise_stim_number
            stimGroup.attrs['song_wavfile'] = repPres.stim_filename
            stimGroup.attrs['noise_wavfile'] = masker.repeated.stim_filename
            stimGroup.attrs['weights'] = 1.0
            stimGroup.attrs['song_start_time'] = float(single.start - masker.start)

            spikeTrials = [masker.spiketimes]

            for trialNum,spikes in enumerate(spikeTrials):
                trialResponseName = '%d' % int(trialNum+1)
                #print '\tCreating trial response with %d spikes: %s' % (len(spikes), trialResponseName)
                trialResponseGroup = stimGroup.create_group(trialResponseName)
                trialResponseGroup.attrs['bg_rate'] = 0.0
                trialResponseGroup.attrs['peri_rate'] = len(spikes) / float(masker.duration)
                if len(spikes) > 0:
                    trialResponseGroup['spike_times'] = np.array(spikes)
                else:
                    trialResponseGroup['spike_times'] = -999



def add_signal_to_noise(signal, noise, snr, start):
    std_signal = signal.std()
    signal = zero_pad(signal, len(noise), start)
    return signal * ((noise.std() * 10 ** (snr / 20)) / std_signal) + noise

def zero_pad(signal, duration, start):
    p_signal = np.zeros(duration)
    p_signal[start: start + len(signal)] = signal
    return p_signal

def get_masked_stims(masked_resps):

    snr = 3.0
    masked_stims = []
    for ii, resp in enumerate(masked_resps):
        masked_stims.append([])
        signal = resp.presentation.stim.file.get_data()
        sr = resp.presentation.stim.file.get_samplerate()
        for sresp in resp.singles:
            if sresp.contained() is None:
                masked_stims.pop()
                break
            nresp = sresp.contained()
            start = sresp.start - nresp.start
            noise = nresp.repeated.presentation.stim.file.get_data()
            start = np.floor(float(start) * sr)
            masked_stims[ii].append(add_signal_to_noise(signal, noise, snr, start))

    return masked_stims

def spectrogram(sound, sample_rate=16000, increment=0, window_length=0, n_std=6, f_band=125, spec_sample_rate=1000):

    desired_sample_rate = np.ceil(np.double(sample_rate) / spec_sample_rate) * spec_sample_rate

    if desired_sample_rate > sample_rate:
        sound = signal.resample(sound, int(float(len(sound)) / sample_rate * desired_sample_rate))
        sample_rate = desired_sample_rate

    if window_length == 0:
        twindow = n_std / (f_band * 2.0 * np.pi)
        window_length = np.fix(twindow * sample_rate)
        window_length = np.fix(window_length / 2.0) * 2.0

    if increment == 0:
        increment = np.floor(sample_rate / spec_sample_rate)

    #Ensure window has an even length
    if np.mod(window_length, 2) == 1:
        window_length = window_length + 1

    #Pad the sound with zeros
    p_sound = np.zeros((len(sound) + window_length,))
    p_sound[window_length / 2: window_length / 2 + len(sound)] = sound
    p_length = len(p_sound)

    frame_count = int(np.floor((p_length - window_length) / increment))

    #Generate Gaussian window
    n_std = 6
    window = np.exp(-.5 * ((np.arange(window_length) - (window_length + 1) / 2.) ** 2 / (window_length / n_std) ** 2))

    fft_sound = np.zeros((window_length / 2 + 1, frame_count))

    for i in range(frame_count):
        start = i * increment
        last = start + window_length
        f = np.zeros((window_length,))
        f = window * p_sound[start: last]
        spec_slice = np.fft.fft(f)
        fft_sound[:, i] = abs(spec_slice[: window_length / 2 + 1])


    freq = abs(np.fft.fftfreq(int(window_length))[: window_length / 2 + 1]) * sample_rate
    time = np.arange(fft_sound.shape[1]) * increment / sample_rate

    return fft_sound, time, freq


def log_spec(spec, dBNoise=80):

    spec = 20 * np.log10(spec / spec.max())
    spec = np.maximum(spec, -dBNoise)
    return spec

def normalize_spec(spec):

    return spec

def conv_activations_2d(stims, filts):

    act = list()
    if type(stims) is not list:
        stims = [stims]

    if type(filts) is not list:
        filts = [filts]

    if len(stims[0].shape) == 1:
        print "Use conv_activations_1d()!"
        return

    for ii, stim in enumerate(stims):
        stim_act = np.zeros((len(filts), stim.shape[1] - 1))
        for jj, filt in enumerate(filts):
            filt_act = np.zeros(stim.shape[1] - 1)
            for kk in range(filt.shape[1]):
                at = np.dot(stim.T, filt[:,kk])
                filt_act[kk:] = filt_act[kk:] + at[:len(filt_act) - kk]
            stim_act[jj] = filt_act
        act.append(stim_act)

    return act