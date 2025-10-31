(async () => {
  const sampleRate = 44100;

  // voice cloning and gender variables
  let genderShift = 0; // -100 (male) to 100 (female)
  let clonedRatios = null; // ratios for formant scaling from cloned voice
  let clonedTimbre = { harmonics: 27, duty: 0.55 }; // default
  const standardFormants = [700, 1220, 2600]; // standard formants for 'a'

  // adjust formants for gender and cloning
  const getAdjustedFormants = (baseF) => {
    let f = [...baseF];
    if (clonedRatios) {
      // apply cloned voice ratios
      f = f.map((freq, i) => freq * clonedRatios[i]);
    }
    // gender shift: scale frequencies (female higher, male lower)
    const scale = 1 + (genderShift / 100) * 0.3; // Â±30%
    f = f.map(freq => freq * scale);
    return f;
  };

  // extract voice parameters from uploaded audio
  const extractVoiceParameters = (audioBuffer) => {
    // Improved extraction: analyze spectrum with time-averaging for better formant detection
    const sampleRate = audioBuffer.sampleRate;
    const data = audioBuffer.getChannelData(0);
    const fftSize = 2048;
    const analyser = new AnalyserNode(new (window.AudioContext || window.webkitAudioContext)(), { fftSize });
    const bufferLength = analyser.frequencyBinCount;
    const dataArray = new Uint8Array(bufferLength);

    // Time-average the spectrum over the entire buffer for sustained vowel assumption
    const avgSpectrum = new Uint8Array(bufferLength);
    const windowSize = Math.floor(data.length / 10); // divide into 10 windows
    for (let w = 0; w < 10; w++) {
      const start = w * windowSize;
      const end = Math.min((w + 1) * windowSize, data.length);
      analyser.getByteFrequencyData(dataArray);
      for (let i = 0; i < bufferLength; i++) {
        avgSpectrum[i] += dataArray[i];
      }
    }
    for (let i = 0; i < bufferLength; i++) {
      avgSpectrum[i] /= 10;
    }

    // Improved peak finding with parabolic interpolation for sub-bin accuracy
    const findPeak = (startHz, endHz) => {
      const startBin = Math.floor(startHz / (sampleRate / 2) * bufferLength);
      const endBin = Math.floor(endHz / (sampleRate / 2) * bufferLength);
      let maxVal = 0;
      let maxBin = startBin;
      for (let i = startBin; i <= endBin && i < bufferLength; i++) {
        if (avgSpectrum[i] > maxVal) {
          maxVal = avgSpectrum[i];
          maxBin = i;
        }
      }
      // Parabolic interpolation for better accuracy
      if (maxBin > startBin && maxBin < endBin) {
        const a = avgSpectrum[maxBin - 1];
        const b = avgSpectrum[maxBin];
        const c = avgSpectrum[maxBin + 1];
        const p = 0.5 * (a - c) / (a - 2 * b + c);
        maxBin += p;
      }
      return maxBin / bufferLength * (sampleRate / 2);
    };

    const f1 = findPeak(400, 1000);
    const f2 = findPeak(1000, 2000);
    const f3 = findPeak(2000, 3500);

    const extractedFormants = [f1 || 700, f2 || 1220, f3 || 2600];
    clonedRatios = extractedFormants.map((f, i) => f / standardFormants[i]);
    clonedTimbre = { harmonics: 27, duty: 0.55 }; // keep default for now
  };

  // phoneme -> formant settings
  const phonemeMap = {
    a: { f: [700, 1220, 2600], voiced: true },
    aa: { f: [900, 1300, 2650], voiced: true },
    e: { f: [500, 2300, 3000], voiced: true },
    i: { f: [350, 2200, 2900], voiced: true },
    ee: { f: [285, 2275, 2900], voiced: true },
    I: { f: [440, 1200, 2700], voiced: true },
    o: { f: [400, 1000, 2600], voiced: true },
    u: { f: [325, 700, 2530], voiced: true },
    y: { f: [300, 2000, 2800], voiced: true },
    w: { f: [400, 1000, 2200], voiced: true },
    r: { f: [450, 1300, 1700], voiced: true },
    l: { f: [500, 1400, 2600], voiced: true },

    // consonants & sibilants (some flagged voiced/unvoiced)
    h: { f: [800, 1800, 3200], breathy: true, amp: 0.9, voiced: false },
    s: { f: [3000, 5000, 7000], breathy: true, amp: 1.4, voiced: false },
    z: { f: [3000, 4500, 6000], breathy: true, amp: 1.4, voiced: true },
    t: { f: [1000, 4500, 7000], burst: true, amp: 0.6, voiced: false },
    d: { f: [600, 4000, 6500], burst: true, amp: 0.6, voiced: false, short: true },
    k: { f: [1200, 2000, 3200], burst: true, short: true, voiced: false },
    g: { f: [900, 1700, 2700], burst: true, voiced: true, short: true },
    // nasals: flagged nasal: true
    n: { f: [300, 1300, 2500], voiced: true, nasal: true },
    m: { f: [250, 1100, 2100], voiced: true, nasal: true },
    b: { f: [800, 1500, 2400], burst: true, voiced: false, short: true },
    p: { f: [1000, 1800, 2700], burst: true, short: true, voiced: false },
    f: { f: [1200, 3000, 5000], breathy: true, voiced: false },
    v: { f: [1200, 2800, 4800], breathy: true, voiced: true },
    th: { f: [1200, 2200, 3500], breathy: true, voiced: false },
    sh: { f: [2500, 3500, 5000], breathy: true, voiced: false },
    ch: { f: [2000, 3000, 4500], breathy: true, burst: true, voiced: false },
    uh: { f: [690, 995, 2600], voiced: true },

    // added phones / fallbacks
    er: { f: [600, 1250, 2500], voiced: true },
    j: { f: [500, 1600, 2600], voiced: true },
    ng: { f: [250, 900, 2000], voiced: true, nasal: true },
    oy: { f: [450, 900, 2600], voiced: true },
    aw: { f: [700, 1300, 2600], voiced: true },
    rest: { f: [0, 0, 0], voiced: false, burst: false, short: false, breathy: false }
  };

  // --- CMUdict-based G2P + ARPAbet -> phonemeKey converter ---
  let CMUDICT = null;
  const CMUDICT_URL = "https://raw.githubusercontent.com/cmusphinx/cmudict/master/cmudict.dict";
  let useCMUDict = false;

  async function loadCMUDict(url = CMUDICT_URL) {
    if (CMUDICT) return CMUDICT;
    const res = await fetch(url);
    if (!res.ok) throw new Error("failed to fetch CMUdict: " + res.status);
    const txt = await res.text();
    const map = new Map();
    const lines = txt.split(/\r?\n/);
    for (const line of lines) {
      if (!line || line.startsWith(";;;")) continue;
      const m = line.match(/^([A-Z'()]+)\s+(.*)$/);
      if (!m) continue;
      const wordRaw = m[1].replace(/\(\d+\)$/, "");
      const word = wordRaw.toLowerCase();
      const arp = m[2].trim().split(/\s+/);
      if (!map.has(word)) map.set(word, arp);
    }
    CMUDICT = map;
    return CMUDICT;
  }

  const arpabetMap = {
    "AA":["aa"],"AE":["a"],"AH":["uh"],"AO":["o"],"AW":["aw"],
    "AY":["aa"],"EH":["e"],"ER":["er"],"EY":["ee"],"IH":["I"],
    "IY":["ee"],"OW":["o"],"OY":["oy"],"UH":["uh"],"UW":["u"],
    "B":["b"],"CH":["ch"],"D":["d"],"DH":["th"],"F":["f"],"G":["g"],
    "HH":["h"],"JH":["j"],"K":["k"],"L":["l"],"M":["m"],"N":["n"],
    "NG":["ng"],"P":["p"],"R":["r"],"S":["s"],"SH":["sh"],"T":["t"],
    "TH":["th"],"V":["v"],"W":["w"],"Y":["y"],"Z":["z"],"ZH":["z"]
  };

  function arpabetToKeys(arpArr) {
    const out = [];
    for (let token of arpArr) {
      token = token.replace(/\d+$/, "");
      const mapped = arpabetMap[token];
      if (mapped) out.push(...mapped);
      else out.push(token.toLowerCase());
    }
    return out;
  }

  async function textToPhonemes(word) {
    if (!word || !word.trim()) return [];
    const w = word.toLowerCase().replace(/[^a-z']/g, "");
    if (w.length === 0) return [];

    if (useCMUDict) {
      try {
        await loadCMUDict();
        if (CMUDICT && CMUDICT.has(w)) {
          const arp = CMUDICT.get(w).slice();
          return arpabetToKeys(arp);
        }
      } catch (e) {
        console.warn("CMUdict load/lookup failed:", e);
      }
    }

    // Greedy fallback (multigraphs + contextual c/g)
    const multigraphRules = [
      ["tion", ["sh","uh","n"]], ["tch", ["ch"]], ["ch", ["ch"]], ["sh", ["sh"]],
      ["ph", ["f"]], ["ng", ["ng"]], ["qu", ["k","w"]], ["wh", ["w"]],
      ["kn", ["n"]], ["wr", ["r"]], ["ee", ["ee"]], ["ea", ["ee"]],
      ["ai", ["aa"]], ["ay", ["aa"]], ["oa", ["o"]], ["oo", ["u"]],
      ["ow", ["o"]], ["oy", ["oy"]], ["au", ["aw"]], ["ough", ["o"]]
    ];
    multigraphRules.sort((a,b) => b[0].length - a[0].length);

    const single = {
      a:"a", b:"b", c:"c", d:"d", e:"e", f:"f", g:"g", h:"h",
      i:"i", j:"j", k:"k", l:"l", m:"m", n:"n", o:"o", p:"p",
      q:"q", r:"r", s:"s", t:"t", u:"u", v:"v", w:"w", x:"x", y:"y", z:"z"
    };

    const out = [];
    let i = 0;
    while (i < w.length) {
      let matched = false;
      for (const [pat, phon] of multigraphRules) {
        if (w.startsWith(pat, i)) {
          out.push(...phon);
          i += pat.length;
          matched = true;
          break;
        }
      }
      if (matched) continue;
      if (w[i] === "'") { i++; continue; }
      const ch = w[i];
      if (ch === "c") {
        const next = w[i+1] || "";
        out.push(...("eiy".includes(next) ? ["s"] : ["k"]));
        i++; continue;
      }
      if (ch === "g") {
        const next = w[i+1] || "";
        out.push(...("eiy".includes(next) ? ["j"] : ["g"]));
        i++; continue;
      }
      if (ch === "x") { out.push("k","s"); i++; continue; }
      if (single[ch]) out.push(single[ch]); else out.push(ch);
      i++;
    }
    return out;
  }

  // helper: convert musical note to frequency, returns numeric freq if note parse fails
  const noteToFreq = (note) => {
    const A4 = 440;
    const semitoneMap = {
      C: -9, "C#": -8, Db: -8,
      D: -7, "D#": -6, Eb: -6,
      E: -5, F: -4, "F#": -3, Gb: -3,
      G: -2, "G#": -1, Ab: -1,
      A: 0, "A#": 1, Bb: 1, B: 2,
    };
    const match = /^([A-Ga-g])([#b]?)(\d)$/.exec(note);
    if (!match) return parseFloat(note);
    const [, base, accidental, octave] = match;
    const name = base.toUpperCase() + (accidental || "");
    const semis = semitoneMap[name] ?? 0;
    const semitonesFromA4 = (parseInt(octave) - 4) * 12 + semis;
    return A4 * Math.pow(2, semitonesFromA4 / 12);
  };

  // parseInput async: supports grid types and uses textToPhonemes for words
  const parseInput = async (text, beatLen, gridType = "beats", stepsPerBeat = 4) => {
    const regex = /([a-zA-Z']+)\s*<([\w#b]+)\s*,\s*([\d.]+)>/gi;
    const result = [];
    let match;
    while ((match = regex.exec(text)) !== null) {
      const tokenRaw = match[1];
      const token = tokenRaw.toLowerCase();
      const pitchRaw = match[2];
      const units = parseFloat(match[3]);

      // compute duration from units based on grid type
      let dur = 0;
      if (gridType === "beats") dur = units * beatLen;
      else if (gridType === "steps") dur = (units / stepsPerBeat) * beatLen;
      else if (gridType === "seconds") dur = units;
      else dur = units * beatLen;

      // exact phoneme key match
      if (phonemeMap[token]) {
        let p = { key: token, ...phonemeMap[token], d: dur, pitch: noteToFreq(pitchRaw) };
        p.f = getAdjustedFormants(p.f);
        result.push(p);
        continue;
      }

      // otherwise treat token as a word and convert to phoneme keys
      const parts = await textToPhonemes(tokenRaw);
      if (!parts || parts.length === 0) continue;
      const partDur = dur / Math.max(1, parts.length);
      for (const part of parts) {
        const key = part.toLowerCase();
        if (phonemeMap[key]) {
          let p = { key, ...phonemeMap[key], d: partDur, pitch: noteToFreq(pitchRaw) };
          p.f = getAdjustedFormants(p.f);
          result.push(p);
        } else {
          const fallbackMap = {
            a: "a", e: "e", i: "i", o: "o", u: "u",
            h: "h", n: "n", m: "m", s: "s", t: "t", d: "d", r: "r", l: "l",
            b: "b", p: "p", k: "k", g: "g", f: "f", v: "v", y: "y", j: "j",
            ng: "ng", er: "er", oy: "oy", aw: "aw"
          };
          const fm = (fallbackMap[key] || "rest");
          let p = { key: fm, ...phonemeMap[fm], d: partDur, pitch: noteToFreq(pitchRaw) };
          p.f = getAdjustedFormants(p.f);
          result.push(p);
        }
      }
    }
    return result;
  };

  // noise buffer for breathy / sibilant sounds
  const createWhisperNoiseBuffer = (ctx, duration) => {
    const len = Math.max(1, Math.floor(duration * ctx.sampleRate));
    const buffer = ctx.createBuffer(1, len, ctx.sampleRate);
    const data = buffer.getChannelData(0);
    for (let i = 0; i < len; i++) data[i] = (Math.random() * 2 - 1) * 0.25;
    return buffer;
  };

  // synthesize with nasal-aware transitions and humanizing features
  const synthesize = async (ctx, phonemeSeq, mode, vibFreq, vibDepth, vibDelay, morphTime = 0.05, morphEnabled = true, slideTime = 0.08, persistentVib = true) => {
    // Preprocess phonemeSeq for humanizing features
    const processedSeq = [];
    for (let i = 0; i < phonemeSeq.length; i++) {
      const p = phonemeSeq[i];
      const prev = processedSeq.length > 0 ? processedSeq[processedSeq.length - 1] : null;
      const hasPrevVoiced = prev && prev.voiced;

      // Feature 3: If b or p and no previous voiced, prepend m
      if ((p.key === 'b' || p.key === 'p') && !hasPrevVoiced) {
        const mDur = Math.min(0.05, p.d * 0.1); // short m
        const mPhoneme = { key: 'm', ...phonemeMap['m'], d: mDur, pitch: p.pitch };
        processedSeq.push(mPhoneme);
        p.d -= mDur; // trim p/b duration
      }

      // Feature 2: If t,b,p,k after s or z, add silence delay
      if ((p.key === 't' || p.key === 'b' || p.key === 'p' || p.key === 'k') && prev && (prev.key === 's' || prev.key === 'z')) {
        const silenceDur = Math.min(0.02, prev.d * 0.1); // slight silence
        prev.d -= silenceDur; // trim s/z
        p.d -= silenceDur; // trim consonant
        // silence is implicit by delaying start
      }

      processedSeq.push(p);
    }

    const totalDuration = processedSeq.reduce((a,p) => a + p.d, 0);
    const endTime = totalDuration;

    // master output
    const master = ctx.createGain();
    master.gain.value = 1;

    // shared vowel voice filters (3 formants) and per-filter gains (so we can dip vowel energy)
    const numFormants = 3;
    const voiceFilters = [];
    const voiceGains = [];
    for (let i = 0; i < numFormants; i++) {
      const f = ctx.createBiquadFilter();
      f.type = "bandpass";
      f.Q.value = 15;
      f.frequency.value = 0;
      const g = ctx.createGain();
      g.gain.value = 1; // normal full vowel energy
      f.connect(g);
      g.connect(master);
      voiceFilters.push(f);
      voiceGains.push(g);
    }

    // Feature 4: Sine wave for bass/timbre backing the pulse train, affected by formant filters
    const bassOsc = ctx.createOscillator();
    bassOsc.type = "sine";
    const bassGain = ctx.createGain();
    bassGain.gain.value = 0.1; // subtle bass
    // Route bass through formant filters for integration
    for (let idx = 0; idx < numFormants; idx++) {
      bassGain.connect(voiceFilters[idx]);
    }
    // Will set frequency and start/stop later

    // shared noise filters for breathy vowel/noise path
    const sharedNoiseFilters = [];
    for (let i = 0; i < numFormants; i++) {
      const f = ctx.createBiquadFilter();
      f.type = "bandpass";
      f.Q.value = 15;
      f.frequency.value = 0;
      f.connect(master);
      sharedNoiseFilters.push(f);
    }

    // persistent vibrato LFO
    let lfo = null;
    let lfoGain = null;
    const vibActive = vibFreq > 0 && vibDepth > 0 && persistentVib && mode !== "whisper";
    if (vibActive) {
      lfo = ctx.createOscillator();
      lfo.type = "sine";
      lfo.frequency.value = vibFreq;
      lfoGain = ctx.createGain();
      lfoGain.gain.value = vibDepth; // Hz depth
      lfo.connect(lfoGain);
      lfo.start(vibDelay);
      lfo.stop(endTime + 1);
    }

    const setFilterNow = (filter, t, val) => {
      try { filter.frequency.setValueAtTime(val, t); } catch (e) {}
    };

    // find next voiced phoneme (skip consonants/unvoiced) and return it
    const findNextVoiced = (i) => {
      for (let j = i + 1; j < phonemeSeq.length; j++) {
        const p = phonemeSeq[j];
        if (p && p.voiced) return p;
      }
      return null;
    };

    // create a short nasal noise burst scheduled in the morph window
    const scheduleNasalBurst = (morphStart, morphEnd, depth = 0.6) => {
      // bandpass centered low (~250-400Hz) for nasal resonance
      const nasalFilter = ctx.createBiquadFilter();
      nasalFilter.type = "bandpass";
      nasalFilter.Q.value = 8;
      nasalFilter.frequency.value = 300; // center freq for general nasal flavour

      const nasalGain = ctx.createGain();
      nasalGain.gain.value = 0;
      nasalFilter.connect(nasalGain).connect(master);

      const dur = Math.max(0.001, morphEnd - morphStart);
      const src = ctx.createBufferSource();
      src.buffer = createWhisperNoiseBuffer(ctx, dur + 0.02);

      // envelope: ramp in at morphStart, ramp out at morphEnd
      nasalGain.gain.setValueAtTime(0, morphStart - 0.001);
      nasalGain.gain.linearRampToValueAtTime(depth, morphStart + Math.min(0.02, dur * 0.25));
      nasalGain.gain.setValueAtTime(depth, morphEnd - Math.min(0.02, dur * 0.25));
      nasalGain.gain.linearRampToValueAtTime(0, morphEnd + 0.01);

      src.connect(nasalFilter);
      src.start(morphStart);
      src.stop(morphEnd + 0.02);
    };

    // consonant noise generator with subtle fade-in
    const playConsonantNoise = (t, d, fArr, amp = 1) => {
      if (d <= 0 || !fArr || fArr.length === 0) return;
      const src = ctx.createBufferSource();
      src.buffer = createWhisperNoiseBuffer(ctx, d);
      const consonantGain = ctx.createGain();
      consonantGain.gain.setValueAtTime(0, t);
      consonantGain.gain.linearRampToValueAtTime(amp, t + 0.01); // subtle fade-in over 0.01s
      const consonantFilters = fArr.map(freq => {
        const bf = ctx.createBiquadFilter();
        bf.type = "bandpass";
        bf.Q.value = 20;
        bf.frequency.value = freq || 0;
        bf.connect(master);
        return bf;
      });
      src.connect(consonantGain);
      for (const bf of consonantFilters) consonantGain.connect(bf);
      src.start(t);
      src.stop(t + d + 0.005);
    };

    // main play routine
    const play = (t, d, f, opt = {}, nextVoiced = null, immediateNext = null) => {
      const { breathy = false, amp = 1, burst = false, voiced = true, pitch = 220, short = false } = opt;
      if (d <= 0) return;

      // If current phoneme is consonant/unvoiced or purely breathy/burst -> produce consonant-filtered noise
      if (!voiced || opt.breathy || opt.burst) {
        playConsonantNoise(t, d, f, amp);
        // If it's not voiced at all, we're done
        if (!voiced) return;
      }

      const shouldVoice = voiced && mode !== "whisper";
      if (!shouldVoice) return;

      // determine morph behaviour: skip consonants for target as before
      const hasVoicedTarget = morphEnabled && morphTime > 0 && nextVoiced && nextVoiced.voiced;

      // if target is nasal, schedule a nasal transition instead of morphing formants toward nasal targets
      const targetIsNasal = hasVoicedTarget && !!nextVoiced.nasal;
      const morphStart = t + Math.max(0, d - morphTime);
      const morphEnd = t + d;

      if (targetIsNasal) {
        // keep vowel formants at current values, but dip vowel gain slightly and add a nasal burst
        for (let idx = 0; idx < numFormants; idx++) {
          const curVal = (f && f[idx]) || 0;
          // hold formant frequency across the phoneme (no ramp to nasal targets)
          setFilterNow(voiceFilters[idx], t, curVal);
          setFilterNow(sharedNoiseFilters[idx], t, curVal);

          // schedule vowel gain dip (small reduction)
          try {
            voiceGains[idx].gain.cancelScheduledValues(t);
            voiceGains[idx].gain.setValueAtTime(1, t);
            // dip to e.g. 0.65 during morph window, then back to 1
            voiceGains[idx].gain.setValueAtTime(1, morphStart - 0.001);
            voiceGains[idx].gain.linearRampToValueAtTime(0.65, morphStart + Math.min(0.02, morphTime * 0.5));
            voiceGains[idx].gain.setValueAtTime(0.65, morphEnd - Math.min(0.02, morphTime * 0.5));
            voiceGains[idx].gain.linearRampToValueAtTime(1, morphEnd + 0.01);
          } catch (e) {}
        }
        // schedule nasal burst (filtered noise) in morph window to add "nah"/nasal timbre
        scheduleNasalBurst(morphStart, morphEnd, Math.max(0.2, 0.6 * (opt.amp || 1)));
      } else if (hasVoicedTarget) {
        // normal morph-to-next voiced formants
        for (let idx = 0; idx < numFormants; idx++) {
          const curVal = (f && f[idx]) || 0;
          const nextVal = (nextVoiced.f && nextVoiced.f[idx]) || 0;
          try {
            voiceFilters[idx].frequency.setValueAtTime(curVal, t);
            voiceFilters[idx].frequency.setValueAtTime(curVal, morphStart);
            voiceFilters[idx].frequency.linearRampToValueAtTime(nextVal, morphEnd);
          } catch (e) {}
          try {
            sharedNoiseFilters[idx].frequency.setValueAtTime(curVal, t);
            sharedNoiseFilters[idx].frequency.setValueAtTime(curVal, morphStart);
            sharedNoiseFilters[idx].frequency.linearRampToValueAtTime(nextVal, morphEnd);
          } catch (e) {}
        }
      } else {
        // no voiced target (or morph disabled) -> set filters to current f immediately
        for (let idx = 0; idx < numFormants; idx++) {
          const val = (f && f[idx]) || 0;
          setFilterNow(voiceFilters[idx], t, val);
          setFilterNow(sharedNoiseFilters[idx], t, val);
          // ensure gains are at full value
          try { voiceGains[idx].gain.setValueAtTime(1, t); } catch (e) {}
        }
      }

      // oscillator (voiced) - custom glottal source approximating Liljencrants-Fant model
      const osc = ctx.createOscillator();

      const numHarmonics = 27;
      const real = new Float32Array(numHarmonics);
      const imag = new Float32Array(numHarmonics);

      // Approximate LF glottal flow derivative spectrum: amplitudes decrease as 1/n^2 for more natural harmonics
      for (let n = 1; n < numHarmonics; n++) {
        real[n] = 1 / (n * n);
      }

      const glottalWave = ctx.createPeriodicWave(real, imag);
      osc.setPeriodicWave(glottalWave);
      
      const pitchParam = osc.frequency;
      pitchParam.setValueAtTime(pitch, t);

      // vibrato
      if (vibActive && lfoGain) {
        lfoGain.connect(pitchParam);
      } else if (!persistentVib && vibFreq > 0 && vibDepth > 0) {
        const lfol = ctx.createOscillator();
        const lfoGl = ctx.createGain();
        lfol.type = "sine";
        lfol.frequency.setValueAtTime(vibFreq, t);
        lfoGl.gain.setValueAtTime(vibDepth, t);
        lfol.connect(lfoGl).connect(pitchParam);
        lfol.start(t + vibDelay);
        lfol.stop(t + d + 0.02);
      }

      // slide: only when immediate next (direct next) is voiced and slideTime > 0
      const canSlide = slideTime > 0 && immediateNext && immediateNext.voiced && immediateNext.pitch && voiced;
      if (canSlide) {
        const rampStart = Math.max(t, t + d - slideTime);
        pitchParam.setValueAtTime(pitch, t);
        pitchParam.setValueAtTime(pitch, rampStart);
        pitchParam.linearRampToValueAtTime(immediateNext.pitch, t + d);
      } else {
        pitchParam.setValueAtTime(pitch, t);
      }

      // route oscillator through per-formant filters/gains
      const oscGain = ctx.createGain();
      const fadeTime = 0.01;
      oscGain.gain.setValueAtTime(0, t);
      oscGain.gain.linearRampToValueAtTime(1 * amp, t + fadeTime);
      oscGain.gain.setValueAtTime(1 * amp, t + d - fadeTime);
      oscGain.gain.linearRampToValueAtTime(0, t + d);
      for (let idx = 0; idx < numFormants; idx++) {
        osc.connect(oscGain).connect(voiceFilters[idx]);
      }

      // If breathy component, route a noise source via sharedNoiseFilters as well
      if (opt.breathy) {
        const noiseSrc = ctx.createBufferSource();
        noiseSrc.buffer = createWhisperNoiseBuffer(ctx, d);
        for (const nf of sharedNoiseFilters) noiseSrc.connect(nf);
        noiseSrc.start(t);
        noiseSrc.stop(t + d + 0.005);
      }

      osc.start(t);
      osc.stop(t + d + fadeTime);
    };

    // initialize voice filters and bass to first voiced phoneme
    if (processedSeq.length > 0) {
      const firstVoiced = processedSeq.find(p => p.voiced) || processedSeq[0];
      if (firstVoiced && firstVoiced.f) {
        for (let idx = 0; idx < numFormants; idx++) {
          const v = firstVoiced.f[idx] || 0;
          setFilterNow(voiceFilters[idx], 0, v);
          setFilterNow(sharedNoiseFilters[idx], 0, v);
        }
        // set bass to octave below
        bassOsc.frequency.setValueAtTime(firstVoiced.pitch / 2, 0);
      }
    }
    // start bass osc
    bassOsc.start(0);
    bassOsc.stop(endTime);

    // continuous oscillator for voiced sequences
    let currentOsc = null;
    let oscGain = null;
    let lastVoicedEnd = 0;
    let currentAmp = 1;

    // schedule phonemes
    let t = 0;
    for (let i = 0; i < processedSeq.length; i++) {
      const p = processedSeq[i];
      const immediateNext = processedSeq[i + 1] || null;
      const nextVoiced = findNextVoiced(i);

      if (p.voiced && mode !== "whisper") {
        if (!currentOsc) {
          // start new oscillator for voiced sequence
          currentOsc = ctx.createOscillator();
          const numHarmonics = 27;
          const real = new Float32Array(numHarmonics);
          const imag = new Float32Array(numHarmonics);
          for (let n = 1; n < numHarmonics; n++) {
            real[n] = 1 / (n * n);
          }
          const glottalWave = ctx.createPeriodicWave(real, imag);
          currentOsc.setPeriodicWave(glottalWave);

          oscGain = ctx.createGain();
          const fadeTime = 0.01;
          const amp = isNaN(p.amp) ? 1 : p.amp;
          oscGain.gain.setValueAtTime(0, t);
          oscGain.gain.linearRampToValueAtTime(2 * amp, t + fadeTime);
          currentAmp = amp;
          for (let idx = 0; idx < numFormants; idx++) {
            currentOsc.connect(oscGain).connect(voiceFilters[idx]);
          }
          currentOsc.start(t);
        }

        // schedule pitch
        const pitchParam = currentOsc.frequency;
        pitchParam.setValueAtTime(p.pitch, t);

        // vibrato
        if (vibActive && lfoGain) {
          lfoGain.connect(pitchParam);
        } else if (!persistentVib && vibFreq > 0 && vibDepth > 0) {
          const lfol = ctx.createOscillator();
          const lfoGl = ctx.createGain();
          lfol.type = "sine";
          lfol.frequency.setValueAtTime(vibFreq, t);
          lfoGl.gain.setValueAtTime(vibDepth, t);
          lfol.connect(lfoGl).connect(pitchParam);
          lfol.start(t + vibDelay);
          lfol.stop(t + p.d + 0.02);
        }

        // slide
        const canSlide = slideTime > 0 && immediateNext && immediateNext.voiced && immediateNext.pitch && p.voiced;
        if (canSlide) {
          const rampStart = Math.max(t, t + p.d - slideTime);
          pitchParam.setValueAtTime(p.pitch, t);
          pitchParam.setValueAtTime(p.pitch, rampStart);
          pitchParam.linearRampToValueAtTime(immediateNext.pitch, t + p.d);
        } else {
          pitchParam.setValueAtTime(p.pitch, t);
        }

        // update lastVoicedEnd
        lastVoicedEnd = t + p.d;

        // morphing logic
        const hasVoicedTarget = morphEnabled && morphTime > 0 && nextVoiced && nextVoiced.voiced;
        const targetIsNasal = hasVoicedTarget && !!nextVoiced.nasal;
        const morphStart = t + Math.max(0, p.d - morphTime);
        const morphEnd = t + p.d;

        if (targetIsNasal) {
          for (let idx = 0; idx < numFormants; idx++) {
            const curVal = (p.f && p.f[idx]) || 0;
            setFilterNow(voiceFilters[idx], t, curVal);
            setFilterNow(sharedNoiseFilters[idx], t, curVal);
            try {
              voiceGains[idx].gain.cancelScheduledValues(t);
              voiceGains[idx].gain.setValueAtTime(1, t);
              voiceGains[idx].gain.setValueAtTime(1, morphStart - 0.001);
              voiceGains[idx].gain.linearRampToValueAtTime(0.65, morphStart + Math.min(0.02, morphTime * 0.5));
              voiceGains[idx].gain.setValueAtTime(0.65, morphEnd - Math.min(0.02, morphTime * 0.5));
              voiceGains[idx].gain.linearRampToValueAtTime(1, morphEnd + 0.01);
            } catch (e) {}
          }
          scheduleNasalBurst(morphStart, morphEnd, Math.max(0.2, 0.6 * (p.amp || 1)));
        } else if (hasVoicedTarget) {
          for (let idx = 0; idx < numFormants; idx++) {
            const curVal = (p.f && p.f[idx]) || 0;
            const nextVal = (nextVoiced.f && nextVoiced.f[idx]) || 0;
            try {
              voiceFilters[idx].frequency.setValueAtTime(curVal, t);
              voiceFilters[idx].frequency.setValueAtTime(curVal, morphStart);
              voiceFilters[idx].frequency.linearRampToValueAtTime(nextVal, morphEnd);
            } catch (e) {}
            try {
              sharedNoiseFilters[idx].frequency.setValueAtTime(curVal, t);
              sharedNoiseFilters[idx].frequency.setValueAtTime(curVal, morphStart);
              sharedNoiseFilters[idx].frequency.linearRampToValueAtTime(nextVal, morphEnd);
            } catch (e) {}
          }
        } else {
          for (let idx = 0; idx < numFormants; idx++) {
            const val = (p.f && p.f[idx]) || 0;
            setFilterNow(voiceFilters[idx], t, val);
            setFilterNow(sharedNoiseFilters[idx], t, val);
            try { voiceGains[idx].gain.setValueAtTime(1, t); } catch (e) {}
          }
        }

        // breathy
        if (p.breathy) {
          const noiseSrc = ctx.createBufferSource();
          noiseSrc.buffer = createWhisperNoiseBuffer(ctx, p.d);
          for (const nf of sharedNoiseFilters) noiseSrc.connect(nf);
          noiseSrc.start(t);
          noiseSrc.stop(t + p.d + 0.005);
        }

      } else {
        // non-voiced or whisper
        if (currentOsc) {
          // stop the oscillator with fade-out
          const fadeTime = 0.01;
          oscGain.gain.setValueAtTime(2 * currentAmp, lastVoicedEnd - fadeTime);
          oscGain.gain.linearRampToValueAtTime(0, lastVoicedEnd);
          currentOsc.stop(lastVoicedEnd);
          currentOsc = null;
          oscGain = null;
        }

        // play consonant noise
        if (!p.voiced || mode === "whisper") {
          playConsonantNoise(t, p.d, p.f, p.amp);
        }
      }

      t += p.d;
    }

    // stop any remaining oscillator
    if (currentOsc) {
      const fadeTime = 0.01;
      oscGain.gain.setValueAtTime(2 * currentAmp, lastVoicedEnd - fadeTime);
      oscGain.gain.linearRampToValueAtTime(0, lastVoicedEnd);
      currentOsc.stop(lastVoicedEnd);
    }

    master.connect(ctx.destination);
    return await ctx.startRendering();
  };

  // UI: minimal container (controls same as previous iterations)
  document.querySelectorAll(".voice-ui").forEach(e => e.remove());
  const container = document.createElement("div");
  container.className = "voice-ui";
  container.style = "position:fixed;top:20px;left:20px;padding:12px;background:#fff;border:1px solid #ccc;font-family:monospace;z-index:9999;max-width:720px;";
  container.innerHTML = `
    <h3 style="margin:0 0 8px 0">HOSTERS FR SYNTHESIZER</h3>
    <div style="font-size:12px;margin-bottom:8px">HOSTERS VERY HUMAN SYNTHESIZER</div>
    <label style="display:block">Melody / Text:<br/>
      <textarea id="phonemeInput" rows="4" style="width:100%;font-family:monospace;">a <C4,0.4> n <C4,0.2> a <C4,0.4></textarea>
    </label>
    <div style="margin-top:6px">
      <label>Voice mode:
        <select id="voiceMode"><option value="voiced" selected>Voiced</option><option value="whisper">Whisper</option></select>
      </label>
      <label style="margin-left:8px">BPM: <input type="number" id="bpm" value="120" style="width:80px"/></label>
      <label style="margin-left:8px">Grid:
        <select id="gridType"><option value="beats" selected>Beats</option><option value="steps">Steps</option><option value="seconds">Seconds</option></select>
      </label>
      <label id="stepsPerBeatLabel" style="margin-left:8px">Steps/Beat: <input type="number" id="stepsPerBeat" value="4" min="1" style="width:60px"/></label>
    </div>
    <div style="margin-top:6px">
      <label>Vibrato Freq (Hz): <input type="number" id="vibFreq" value="6" step="0.1" style="width:80px"/></label>
      <label style="margin-left:8px">Depth (Hz): <input type="number" id="vibDepth" value="5" step="0.1" style="width:80px"/></label>
      <label style="margin-left:8px">Delay (s): <input type="number" id="vibDelay" value="0.1" step="0.01" style="width:80px"/></label>
      <label style="margin-left:8px"><input type="checkbox" id="persistentVib" checked/> Persistent Vibrato</label>
    </div>
    <div style="margin-top:6px">
      <label>Formant Morph time (s): <input type="number" id="morphTime" value="0.06" step="0.01" style="width:80px"/></label>
      <label style="margin-left:8px"><input type="checkbox" id="enableMorph" checked/> Enable Morph</label>
      <label style="margin-left:8px">Slide time (s): <input type="number" id="slideTime" value="0.08" step="0.01" style="width:80px"/></label>
    </div>
    <div style="margin-top:6px">
      <label><input type="checkbox" id="useCMU"/> Use CMUDict</label>
      <button id="loadCMU" style="margin-left:8px;padding:4px 8px;">Load CMUDict</button>
    </div>
    <div style="margin-top:6px">
      <label>Gender Shift: <input type="range" id="genderShift" min="-100" max="100" value="0" style="width:200px"/></label>
      <span id="genderValue">0</span>
    </div>
    <div style="margin-top:6px">
      <label>Voice Clone: <input type="file" id="voiceFile" accept="audio/*"/></label>
      <div id="cloneControls" style="margin-top:4px;"></div>
    </div>
    <div style="margin-top:8px">
      <button id="synthBtn" style="padding:6px 12px">ðŸŽ¤ Synthesize</button>
      <span id="status" style="margin-left:10px;font-size:12px;color:#444"></span>
    </div>
    <div id="outputControls" style="margin-top:8px"></div>
  `;
  document.body.appendChild(container);

  const gridTypeEl = container.querySelector("#gridType");
  const stepsLabel = container.querySelector("#stepsPerBeatLabel");
  gridTypeEl.addEventListener("change", () => {
    stepsLabel.style.display = gridTypeEl.value === "steps" ? "inline-block" : "none";
  });
  stepsLabel.style.display = gridTypeEl.value === "steps" ? "inline-block" : "none";

  // Gender shift slider event listener
  const genderShiftEl = container.querySelector("#genderShift");
  const genderValueEl = container.querySelector("#genderValue");
  genderShiftEl.addEventListener("input", (e) => {
    genderShift = parseInt(e.target.value);
    genderValueEl.textContent = genderShift;
  });
  genderValueEl.textContent = genderShift; // initial display

  // Voice clone file upload event listener
  const voiceFileEl = container.querySelector("#voiceFile");
  const cloneControlsEl = container.querySelector("#cloneControls");
  voiceFileEl.addEventListener("change", async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    cloneControlsEl.innerHTML = "Extracting voice parameters...";
    try {
      const arrayBuffer = await file.arrayBuffer();
      const audioCtx = new (window.AudioContext || window.webkitAudioContext)();
      const audioBuffer = await audioCtx.decodeAudioData(arrayBuffer);
      extractVoiceParameters(audioBuffer);
      cloneControlsEl.innerHTML = `<button id="cloneBtn">Clone Voice</button> <span>Voice extracted successfully.</span>`;
      const cloneBtn = cloneControlsEl.querySelector("#cloneBtn");
      cloneBtn.addEventListener("click", () => {
        // Toggle cloning: if clonedRatios exists, clear it; else set it
        if (clonedRatios) {
          clonedRatios = null;
          cloneBtn.textContent = "Clone Voice";
          cloneControlsEl.querySelector("span").textContent = "Voice cloning disabled.";
        } else {
          // Re-extract or just enable
          extractVoiceParameters(audioBuffer); // re-extract if needed, but since it's already done, just set
          cloneBtn.textContent = "Disable Clone";
          cloneControlsEl.querySelector("span").textContent = "Voice cloned successfully.";
        }
      });
    } catch (err) {
      cloneControlsEl.innerHTML = "Error extracting voice: " + err.message;
    }
  });

  const synthBtn = container.querySelector("#synthBtn");
  const outputControls = container.querySelector("#outputControls");
  const statusEl = container.querySelector("#status");
  const useCMUChk = container.querySelector("#useCMU");
  const loadCMUBtn = container.querySelector("#loadCMU");

  useCMUChk.addEventListener("change", (e) => useCMUDict = e.target.checked);
  loadCMUBtn.addEventListener("click", async () => {
    try {
      loadCMUBtn.disabled = true; loadCMUBtn.textContent = "Loading...";
      await loadCMUDict();
      loadCMUBtn.textContent = "Loaded âœ“"; useCMUChk.checked = true;
    } catch (err) {
      console.error(err); loadCMUBtn.textContent = "Load failed";
    } finally { loadCMUBtn.disabled = false; }
  });

  synthBtn.onclick = async () => {
    outputControls.innerHTML = ""; statusEl.textContent = "Parsing...";
    try {
      const text = container.querySelector("#phonemeInput").value;
      const mode = container.querySelector("#voiceMode").value;
      const bpm = parseFloat(container.querySelector("#bpm").value) || 120;
      const beatLen = 60 / bpm;
      const gridType = container.querySelector("#gridType").value || "beats";
      const stepsPerBeat = Math.max(1, parseInt(container.querySelector("#stepsPerBeat").value) || 4);
      const vibratoFreq = parseFloat(container.querySelector("#vibFreq").value) || 0;
      const vibratoDepth = parseFloat(container.querySelector("#vibDepth").value) || 0;
      const vibratoDelay = parseFloat(container.querySelector("#vibDelay").value) || 0;
      const persistentVib = !!container.querySelector("#persistentVib").checked;
      const morphTime = Math.max(0, parseFloat(container.querySelector("#morphTime").value) || 0.06);
      const enableMorph = !!container.querySelector("#enableMorph").checked;
      const slideTime = Math.max(0, parseFloat(container.querySelector("#slideTime").value) || 0.08);

      const phonemeSeq = await parseInput(text, beatLen, gridType, stepsPerBeat);
      if (!phonemeSeq || phonemeSeq.length === 0) { statusEl.textContent = "Parsed no phonemes â€” check input."; return; }

      const totalDuration = phonemeSeq.reduce((a,p) => a + p.d, 0);
      statusEl.textContent = `Rendering ${totalDuration.toFixed(2)}s...`;
      const offlineCtx = new OfflineAudioContext(1, Math.ceil(totalDuration * sampleRate) + 128, sampleRate);

      const renderedBuffer = await synthesize(offlineCtx, phonemeSeq, mode, vibratoFreq, vibratoDepth, vibratoDelay, morphTime, enableMorph, slideTime, persistentVib);
      statusEl.textContent = "Done";

      // WAV creation & preview
      const audioData = renderedBuffer.getChannelData(0);
      const wavBuffer = new ArrayBuffer(44 + audioData.length * 2);
      const view = new DataView(wavBuffer);
      const writeString = (offset, str) => {
        for (let i = 0; i < str.length; i++) view.setUint8(offset + i, str.charCodeAt(i));
      };
      writeString(0, "RIFF");
      view.setUint32(4, 36 + audioData.length * 2, true);
      writeString(8, "WAVE");
      writeString(12, "fmt ");
      view.setUint32(16, 16, true);
      view.setUint16(20, 1, true);
      view.setUint16(22, 1, true);
      view.setUint32(24, sampleRate, true);
      view.setUint32(28, sampleRate * 2, true);
      view.setUint16(32, 2, true);
      view.setUint16(34, 16, true);
      writeString(36, "data");
      view.setUint32(40, audioData.length * 2, true);
      for (let i = 0; i < audioData.length; i++) {
        const s = Math.max(-1, Math.min(1, audioData[i]));
        view.setInt16(44 + i * 2, s * 32767, true);
      }
      const blob = new Blob([view], { type: "audio/wav" });
      const url = URL.createObjectURL(blob);
      const dl = document.createElement("a");
      dl.href = url; dl.download = `singer-${Date.now()}.wav`; dl.textContent = "Download WAV";
      dl.style = "display:inline-block;margin-right:8px;margin-top:8px;";
      outputControls.appendChild(dl);

      const playBtn = document.createElement("button");
      playBtn.textContent = "â–¶ Preview"; playBtn.style = "margin-top:8px;padding:6px 12px;";
      playBtn.onclick = () => {
        const actx = new (window.AudioContext || window.webkitAudioContext)();
        const src = actx.createBufferSource();
        src.buffer = renderedBuffer;
        src.connect(actx.destination);
        src.start();
      };
      outputControls.appendChild(playBtn);

    } catch (err) {
      console.error(err);
      statusEl.textContent = "Error: " + (err.message || err);
    }
  };
})();
