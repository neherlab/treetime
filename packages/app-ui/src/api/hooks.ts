import type { ApiClient } from "@neherlab/app-contracts/client";
import {
  queryOptions,
  useMutation,
  useQuery,
  useQueryClient,
  type MutationOptions,
  type QueryClient,
  type UseMutationResult,
  type UseQueryOptions,
  type UseQueryResult,
} from "@tanstack/react-query";

import { useApiContext } from "./context";
import { requestKey, type ApiCallContext } from "./keys";

export type ApiCall<TResult> = (context: ApiCallContext) => Promise<{ data: TResult }>;

export type ApiMutationCall<TVariables, TResult> = (
  context: ApiCallContext,
  variables: TVariables,
) => Promise<{ data: TResult }>;

export type ApiQueryOptions<TResult> = Omit<UseQueryOptions<TResult>, "queryKey" | "queryFn">;

export interface ApiMutationOptions<TVariables, TResult> {
  seed?: (result: TResult, variables: TVariables) => ApiCall<TResult> | undefined;
}

export function useApi<TResult>(
  call: ApiCall<TResult>,
  options: ApiQueryOptions<TResult> = {},
): UseQueryResult<TResult> {
  const { client } = useApiContext();

  return useQuery(apiQueryOptions(client, call, options));
}

export function useApiMutation<TVariables, TResult>(
  call: ApiMutationCall<TVariables, TResult>,
  options: ApiMutationOptions<TVariables, TResult> = {},
): UseMutationResult<TResult, Error, TVariables> {
  const { client } = useApiContext();
  const queryClient = useQueryClient();

  return useMutation(apiMutationOptions(client, queryClient, call, options));
}

export function apiQueryOptions<TResult>(
  client: ApiClient,
  call: ApiCall<TResult>,
  options: ApiQueryOptions<TResult> = {},
) {
  return queryOptions({
    ...options,
    queryKey: requestKey(client, call),
    queryFn: async ({ signal }) => (await call({ client, throwOnError: true, signal })).data,
  });
}

export function apiKey<TResult>(client: ApiClient, call: ApiCall<TResult>) {
  return apiQueryOptions(client, call).queryKey;
}

export function apiMutationOptions<TVariables, TResult>(
  client: ApiClient,
  queryClient: QueryClient,
  call: ApiMutationCall<TVariables, TResult>,
  options: ApiMutationOptions<TVariables, TResult> = {},
): MutationOptions<TResult, Error, TVariables> {
  return {
    mutationFn: async (variables) => (await call({ client, throwOnError: true }, variables)).data,
    onSuccess: (result, variables) => {
      const seed = options.seed?.(result, variables);

      if (seed !== undefined) {
        queryClient.setQueryData(apiKey(client, seed), result);
      }
    },
  };
}
