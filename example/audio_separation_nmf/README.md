# Audio separation using STFT and non-negative matrix factorization

## usage

```console
$ dub run
```

Now, some png and wav are created. You can listen every four notes in test10k.wav are factorized into 0-3.wav.

## results

- mixed waveform

![mixed](mixed.png)

- factorized waveform

![factorized](factorized.png)

- mixed waveform in STFT domain

![mixed-stft](mixed-stft.png)

- factorized waveform and STFT

![factorized-stft](time_basis.png)

## See also

Julia reference implementation https://github.com/r9y9/julia-nmf-ss-toy
