#import urllib.request
import sys
import urllib
# import path
import pkg_resources
import os
from straintypemer import _ROOT
from straintypemer import mlst_urls

def update_mlst_resources():
    resource_path = os.path.join(_ROOT, "mlst_resources/")
    for k in mlst_urls.keys():
        strain_path = os.path.join(resource_path, k)
        if os.path.exists(strain_path) is False:
            os.mkdir(strain_path)

        sys.stderr.write("Updating '{0}' resources\n".format(k))
        for url in mlst_urls[k]:
             sys.stderr.write('\tretrieving: {0}\n'.format(url))
             urllib.urlretrieve(url, os.path.join(strain_path, url.split("/")[-1]  ))

