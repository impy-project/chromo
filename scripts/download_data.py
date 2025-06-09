from chromo.util import _cached_data_dir, get_all_models

urls = set()
for Model in get_all_models():
    if hasattr(Model, "_data_url"):
        url = Model._data_url
        urls.add(url)

for url in urls:
    print(f"Downloading/installing from {url}")  # noqa: T201
    _cached_data_dir(url)
