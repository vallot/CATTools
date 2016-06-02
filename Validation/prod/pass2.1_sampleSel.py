#!/usr/bin/env python

import sys, os
if not os.path.exists("pass2"): os.mkdir("pass2")

from math import sqrt
import json
ds = {}
dsIn = json.loads(open("pass1/dataset.json").read())
for name in dsIn:
    x = dsIn[name]
    title = x["title"]
    safeTitle = title
    for i in " ;:!@#$%^&*()-+=/<>?[]{}|": safeTitle = safeTitle.replace(i,'_')
    if title not in ds:
        ds[title] = {
            'colour':x['colour'],
            'hist':'pass2/central/%s.root' % safeTitle, ## Path to the merged histogram
            'samples':[], ## List of input samples
        }

    ds[title]['samples'].append({
        'type':x['type'],
        'avgWgt':x['avgWgt'], 'nevt':x['nevt'],
        'xsec':x['xsec'], 'lumi':x['lumi'], 'normFactor':x['normFactor'],
        'hist':x['hist'],
    })

open("pass2/dataset.json", "w").write(json.dumps(ds, sort_keys=True, indent=4))
